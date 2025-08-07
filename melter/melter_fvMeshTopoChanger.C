/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "melter_fvMeshTopoChanger.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H" 
#include "addToRunTimeSelectionTable.H"
#include "polyDistributionMap.H" 
#include "OFstream.H"
#include "fileName.H"
#include "OSspecific.H"
#include "mappedFvPatchBaseBase.H"  // 新增：用于获取耦合映射

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    // 添加这两行 - 这是缺失的关键部分！
    defineTypeNameAndDebug(melter, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, melter, fvMesh);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::fvMeshTopoChangers::melter::buildCouplingPointMapping(MeltingDataOutput& outputData) const
{
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    
    if (couplingPatchID_ >= 0 && couplingPatchID_ < pbm.size())
    {
        const polyPatch& solidToGasPatch = pbm[couplingPatchID_];
        
        // 检查是否为映射边界
        const fvPatch& fvp = mesh().boundary()[couplingPatchID_];
        
        if (isA<mappedFvPatchBaseBase>(fvp))
        {
            const mappedFvPatchBaseBase& mapper = 
                refCast<const mappedFvPatchBaseBase>(fvp);
            
            try
            {
                Info<< nl << "=== 耦合面点映射分析 ===" << nl;
                Info<< "固相patch: " << solidToGasPatch.name() 
                    << " (ID: " << couplingPatchID_ << ")" << nl;

                // 获取邻域patch
                const fvPatch& patchNbr = mapper.nbrFvPatch();
                const polyPatch& polyPatchNbr = patchNbr.patch();
                
                Info<< "气相patch: " << polyPatchNbr.name() << nl;
                Info<< "固相面数: " << solidToGasPatch.size() << nl;
                Info<< "气相面数: " << polyPatchNbr.size() << nl;
                
                // 检查面数量是否匹配
                if (solidToGasPatch.size() != polyPatchNbr.size())
                {
                    WarningInFunction
                        << "固相和气相patch的面数量不匹配: "
                        << "固相面数 = " << solidToGasPatch.size()
                        << ", 气相面数 = " << polyPatchNbr.size() << nl
                        << "跳过点映射建立" << endl;
                    return;
                }
                
                if (solidToGasPatch.size() == 0)
                {
                    Info<< "patch为空，跳过点映射建立" << endl;
                    return;
                }
                
                // 获取网格点坐标
                const pointField& solidMeshPoints = mesh().points();
                const pointField& gasMeshPoints = patchNbr.boundaryMesh().mesh().points();
                
                if (solidMeshPoints.size() == 0 || gasMeshPoints.size() == 0)
                {
                    WarningInFunction
                        << "网格点坐标数组为空，跳过点映射建立" << endl;
                    return;
                }
                
                // 统计变量
                label nMappedPoints = 0;
                label nValidMappedPoints = 0;
                scalar totalDistance = 0.0;
                scalar maxDistance = 0.0;
                label nLargeDistancePoints = 0;
                
                const scalar tolerableDistance = 1e-3; // 放宽容差到1毫米
                
                Info<< "开始基于面对应关系的点映射..." << nl;
                
                // 限制处理的面数量，避免处理过多面时出错
                const label maxFacesToProcess = solidToGasPatch.size();
                Info<< "将处理前" << maxFacesToProcess << "个面进行测试..." << nl;
                
                // 修正：直接使用face中的全局点ID，不需要通过meshPoints转换
                for (label facei = 0; facei < maxFacesToProcess; facei++)
                {
                    // 安全检查：确保面索引在有效范围内
                    if (facei >= solidToGasPatch.size() || facei >= polyPatchNbr.size())
                    {
                        WarningInFunction
                            << "面索引" << facei << "超出范围，跳过" << endl;
                        continue;
                    }
                    
                    const face& solidFace = solidToGasPatch[facei];
                    const face& gasFace = polyPatchNbr[facei];
                    
                    // 检查面点数是否匹配
                    if (solidFace.size() != gasFace.size())
                    {
                        WarningInFunction
                            << "面" << facei << "的点数不匹配: "
                            << "固相=" << solidFace.size() 
                            << ", 气相=" << gasFace.size() << endl;
                        continue;
                    }
                    
                    if (solidFace.size() == 0 || gasFace.size() == 0)
                    {
                        WarningInFunction
                            << "面" << facei << "为空面，跳过" << endl;
                        continue;
                    }
                    
                    // 修正：直接使用face中的全局点ID
                    pointField solidFacePoints(solidFace.size());
                    labelList solidFaceGlobalIds(solidFace.size());
                    
                    bool solidFaceValid = true;
                    forAll(solidFace, pointi)
                    {
                        // 直接使用face中的ID作为全局点ID
                        label solidGlobalId = solidFace[pointi];
                        
                        // 检查全局ID是否在有效范围内
                        if (solidGlobalId < 0 || solidGlobalId >= solidMeshPoints.size())
                        {
                            WarningInFunction
                                << "面" << facei << "固相全局点ID" << solidGlobalId 
                                << "超出范围[0," << solidMeshPoints.size()-1 << "]" << endl;
                            solidFaceValid = false;
                            break;
                        }
                        
                        solidFacePoints[pointi] = solidMeshPoints[solidGlobalId];
                        solidFaceGlobalIds[pointi] = solidGlobalId;
                    }
                    
                    if (!solidFaceValid)
                    {
                        WarningInFunction
                            << "面" << facei << "固相点索引无效，跳过" << endl;
                        continue;
                    }
                    
                    // 修正：直接使用face中的全局点ID
                    pointField gasFacePoints(gasFace.size());
                    labelList gasFaceGlobalIds(gasFace.size());
                    
                    bool gasFaceValid = true;
                    forAll(gasFace, pointi)
                    {
                        // 直接使用face中的ID作为全局点ID
                        label gasGlobalId = gasFace[pointi];
                        
                        // 检查全局ID是否在有效范围内
                        if (gasGlobalId < 0 || gasGlobalId >= gasMeshPoints.size())
                        {
                            WarningInFunction
                                << "面" << facei << "气相全局点ID" << gasGlobalId 
                                << "超出范围[0," << gasMeshPoints.size()-1 << "]" << endl;
                            gasFaceValid = false;
                            break;
                        }
                        
                        gasFacePoints[pointi] = gasMeshPoints[gasGlobalId];
                        gasFaceGlobalIds[pointi] = gasGlobalId;
                    }
                    
                    if (!gasFaceValid)
                    {
                        WarningInFunction
                            << "面" << facei << "气相点索引无效，跳过" << endl;
                        continue;
                    }
                    
                    // 对固相面上的每个点，在气相面上找到最近的点
                    forAll(solidFacePoints, solidPointi)
                    {
                        const point& solidPt = solidFacePoints[solidPointi];
                        const label solidGlobalId = solidFaceGlobalIds[solidPointi];
                        
                        // 在气相面的点中搜索最近的点
                        label bestGasPointi = -1;
                        scalar minDistance = GREAT;
                        
                        forAll(gasFacePoints, gasPointi)
                        {
                            scalar distance = mag(solidPt - gasFacePoints[gasPointi]);
                            if (distance < minDistance)
                            {
                                minDistance = distance;
                                bestGasPointi = gasPointi;
                            }
                        }
                        
                        nMappedPoints++;
                        totalDistance += minDistance;
                        maxDistance = max(maxDistance, minDistance);
                        
                        // 详细输出前几个面的前几个点的映射情况
                        if (facei < 2 && solidPointi < 2)
                        {
                            Info<< "    面" << facei << "固相点" << solidPointi 
                                << "(全局" << solidGlobalId << "): " << solidPt << " -> ";
                            
                            if (bestGasPointi >= 0)
                            {
                                const label gasGlobalId = gasFaceGlobalIds[bestGasPointi];
                                const point& gasPt = gasFacePoints[bestGasPointi];
                                Info<< "气相点" << bestGasPointi 
                                    << "(全局" << gasGlobalId << "): " 
                                    << gasPt << " 距离=" << minDistance << nl;
                            }
                            else
                            {
                                Info<< "未找到匹配点!" << nl;
                            }
                        }
                        
                        // 检查距离是否在容差范围内
                        if (minDistance <= tolerableDistance && bestGasPointi >= 0)
                        {
                            const label gasGlobalId = gasFaceGlobalIds[bestGasPointi];
                            
                            // 检查是否已经被其他固相点映射过
                            bool alreadyMapped = false;
                            forAllConstIter(Map<label>, outputData.solidToGasPointMapping, iter)
                            {
                                if (iter() == gasGlobalId)
                                {
                                    alreadyMapped = true;
                                    break;
                                }
                            }
                            
                            if (!alreadyMapped)
                            {
                                // 建立映射关系
                                outputData.solidToGasPointMapping.insert(solidGlobalId, gasGlobalId);
                                nValidMappedPoints++;
                            }
                        }
                        else
                        {
                            nLargeDistancePoints++;
                            if (nLargeDistancePoints <= 5)
                            {
                                WarningInFunction
                                    << "面" << facei << "固相点" << solidGlobalId 
                                    << "到最近气相点的距离过大: " 
                                    << minDistance << " m" << nl
                                    << "固相点坐标: " << solidPt << nl
                                    << "容差: " << tolerableDistance << " m" << endl;
                            }
                        }
                    }
                }
                
                // 输出映射统计信息
                Info<< nl << "基于面对应关系的点映射统计:" << nl;
                Info<< "  检查的固相点数: " << nMappedPoints << nl;
                Info<< "  成功建立映射的点数: " << nValidMappedPoints << nl;
                Info<< "  距离过大的点数: " << nLargeDistancePoints << nl;
                
                if (nMappedPoints > 0)
                {
                    Info<< "  平均映射距离: " << totalDistance/nMappedPoints << " m" << nl;
                }
                Info<< "  最大映射距离: " << maxDistance << " m" << nl;
                Info<< "  容差: " << tolerableDistance << " m" << nl;
                
                // 检查映射质量
                const scalar mappingRatio = scalar(nValidMappedPoints) / max(nMappedPoints, 1);
                
                if (nValidMappedPoints == 0)
                {
                    WarningInFunction
                        << "未能建立任何点映射关系!" << nl
                        << "这可能表明两个patch的几何形状不匹配或点顺序完全不同" << nl
                        << "将跳过点映射，使用原始点索引" << endl;
                }
                else
                {
                    if (maxDistance < 1e-12)
                    {
                        Info<< "  映射质量: 优秀 (完全重合)" << nl;
                    }
                    else if (maxDistance < 1e-9)
                    {
                        Info<< "  映射质量: 很好 (纳米级精度)" << nl;
                    }
                    else if (maxDistance < 1e-6)
                    {
                        Info<< "  映射质量: 良好 (微米级精度)" << nl;
                    }
                    else if (maxDistance < 1e-3)
                    {
                        Info<< "  映射质量: 可接受 (毫米级精度)" << nl;
                    }
                    else
                    {
                        Info<< "  映射质量: 较差 (距离较大)" << nl;
                    }
                }
                
                Info<< "最终建立了 " << outputData.solidToGasPointMapping.size() 
                    << " 个有效的耦合界面点映射关系" << nl;
                Info<< "映射成功率: " << (mappingRatio * 100) << "%" << nl;
                Info<< "=========================" << nl << endl;
            }
            catch (const std::exception& e)
            {
                WarningInFunction
                    << "无法获取邻域patch信息: " << e.what() << nl
                    << "将跳过点映射建立" << endl;
            }
            catch (...)
            {
                WarningInFunction
                    << "发生未知异常，将跳过点映射建立" << endl;
            }
        }
        else
        {
            Info<< nl << "=== 跳过耦合面点映射 ===" << nl;
            Info<< "固相patch " << solidToGasPatch.name() 
                << " 不是映射边界，边界类型: " << fvp.type() << nl
                << "将使用原始点索引，不进行点映射转换" << nl
                << "===========================" << nl << endl;
        }
    }
    else
    {
        WarningInFunction
            << "无效的solid_to_gas边界ID: " << couplingPatchID_ << nl
            << "边界总数: " << pbm.size() << endl;
    }
}

// void Foam::fvMeshTopoChangers::melter::buildCouplingPointMapping(MeltingDataOutput& outputData) const
// {
//     const polyBoundaryMesh& pbm = mesh().boundaryMesh();
//     if (couplingPatchID_ >= 0 && couplingPatchID_ < pbm.size())
//     {
//         const polyPatch& solidToGasPatch = pbm[couplingPatchID_];     
//         // 检查是否为映射边界
//         const fvPatch& fvp = mesh().boundary()[couplingPatchID_];   
//         if (isA<mappedFvPatchBaseBase>(fvp))
//         {
//             const mappedFvPatchBaseBase& mapper = 
//                 refCast<const mappedFvPatchBaseBase>(fvp);          
//             try
//             {
//                 Info<< nl << "=== 耦合面点映射分析（基于patch点遍历）===" << nl;
//                 Info<< "固相patch: " << solidToGasPatch.name() 
//                     << " (ID: " << couplingPatchID_ << ")" << nl;
//                 // 获取邻域patch
//                 const fvPatch& patchNbr = mapper.nbrFvPatch();
//                 const polyPatch& polyPatchNbr = patchNbr.patch();                
//                 Info<< "气相patch: " << polyPatchNbr.name() << nl;
//                 Info<< "固相面数: " << solidToGasPatch.size() << nl;
//                 Info<< "气相面数: " << polyPatchNbr.size() << nl;                
//                 // 检查面数量是否匹配
//                 if (solidToGasPatch.size() != polyPatchNbr.size())
//                 {
//                     WarningInFunction
//                         << "固相和气相patch的面数量不匹配: "
//                         << "固相面数 = " << solidToGasPatch.size()
//                         << ", 气相面数 = " << polyPatchNbr.size() << nl
//                         << "跳过点映射建立" << endl;
//                     return;
//                 }                
//                 if (solidToGasPatch.size() == 0)
//                 {
//                     Info<< "patch为空，跳过点映射建立" << endl;
//                     return;
//                 }               
//                 // 获取网格点坐标
//                 const pointField& solidMeshPoints = mesh().points();
//                 const pointField& gasMeshPoints = patchNbr.boundaryMesh().mesh().points();                
//                 if (solidMeshPoints.size() == 0 || gasMeshPoints.size() == 0)
//                 {
//                     WarningInFunction
//                         << "网格点坐标数组为空，跳过点映射建立" << endl;
//                     return;
//                 }             
//                 // *** 新方法：基于patch点遍历的映射 ***                
//                 // 1. 获取固相patch的所有点（meshPoints返回patch边界上的点的全局ID）
//                 const labelList& solidPatchPoints = solidToGasPatch.meshPoints();
//                 const labelList& gasPatchPoints = polyPatchNbr.meshPoints();               
//                 Info<< "固相patch点数: " << solidPatchPoints.size() << nl;
//                 Info<< "气相patch点数: " << gasPatchPoints.size() << nl;               
//                 if (solidPatchPoints.size() == 0 || gasPatchPoints.size() == 0)
//                 {
//                     WarningInFunction
//                         << "patch点列表为空，跳过点映射建立" << endl;
//                     return;
//                 }              
//                 // 2. 提取气相patch点的坐标，用于快速搜索
//                 pointField gasPatchPointCoords(gasPatchPoints.size());
//                 forAll(gasPatchPoints, i)
//                 {
//                     label gasGlobalId = gasPatchPoints[i];
//                     if (gasGlobalId >= 0 && gasGlobalId < gasMeshPoints.size())
//                     {
//                         gasPatchPointCoords[i] = gasMeshPoints[gasGlobalId];
//                     }
//                     else
//                     {
//                         WarningInFunction
//                             << "气相patch点" << i << "的全局ID " << gasGlobalId 
//                             << " 超出范围[0," << gasMeshPoints.size()-1 << "]" << endl;
//                         return;
//                     }
//                 }               
//                 // 3. 创建已映射点的标记数组（用于避免重复映射）
//                 PackedBoolList gasPatchPointUsed(gasPatchPoints.size(), false);             
//                 // 统计变量
//                 label nMappedPoints = 0;
//                 label nValidMappedPoints = 0;
//                 scalar totalDistance = 0.0;
//                 scalar maxDistance = 0.0;
//                 label nLargeDistancePoints = 0;             
//                 const scalar tolerableDistance = 1e-3; // 容差：1毫米             
//                 // Info<< "开始基于patch点遍历的映射..." << nl;
//                 // Info<< "容差设置: " << tolerableDistance << " m" << nl;             
//                 // 4. 遍历固相patch的每个点
//                 forAll(solidPatchPoints, solidPointi)
//                 {
//                     label solidGlobalId = solidPatchPoints[solidPointi];                   
//                     // 验证固相点ID有效性
//                     if (solidGlobalId < 0 || solidGlobalId >= solidMeshPoints.size())
//                     {
//                         WarningInFunction
//                             << "固相patch点" << solidPointi << "的全局ID " << solidGlobalId 
//                             << " 超出范围[0," << solidMeshPoints.size()-1 << "]" << endl;
//                         continue;
//                     }                 
//                     const point& solidPt = solidMeshPoints[solidGlobalId];                 
//                     // 在气相patch的所有未使用点中搜索最近的点
//                     label bestGasPointi = -1;
//                     scalar minDistance = GREAT;                  
//                     forAll(gasPatchPointCoords, gasPointi)
//                     {
//                         // 跳过已经被映射的点
//                         if (gasPatchPointUsed[gasPointi])
//                         {
//                             continue;
//                         }                   
//                         scalar distance = mag(solidPt - gasPatchPointCoords[gasPointi]);
//                         if (distance < minDistance)
//                         {
//                             minDistance = distance;
//                             bestGasPointi = gasPointi;
//                         }
//                     }                  
//                     nMappedPoints++;                 
//                     if (bestGasPointi >= 0)
//                     {
//                         totalDistance += minDistance;
//                        maxDistance = max(maxDistance, minDistance);                     
//                         // 详细输出前几个点的映射情况
//                         if (solidPointi < 10)
//                         {
//                             label gasGlobalId = gasPatchPoints[bestGasPointi];
//                             const point& gasPt = gasPatchPointCoords[bestGasPointi];                           
//                             Info<< "  固相patch点" << solidPointi 
//                                 << "(全局" << solidGlobalId << "): " << solidPt << " -> "
//                                 << "气相patch点" << bestGasPointi 
//                                 << "(全局" << gasGlobalId << "): " 
//                                 << gasPt << " 距离=" << minDistance << nl;
//                         }                       
//                         // 检查距离是否在容差范围内
//                         if (minDistance <= tolerableDistance)
//                         {
//                             label gasGlobalId = gasPatchPoints[bestGasPointi];                            
//                             // 建立映射关系
//                             outputData.solidToGasPointMapping.insert(solidGlobalId, gasGlobalId);
//                             nValidMappedPoints++;                           
//                             // *** 关键改进：标记此气相点为已使用，避免重复映射 ***
//                             gasPatchPointUsed[bestGasPointi] = true;                           
//                             // if (debugLevel_ >= 3)
//                             // {
//                             //     Info<< "    建立映射: 固相" << solidGlobalId 
//                             //         << " -> 气相" << gasGlobalId 
//                             //         << " (距离=" << minDistance << ")" << endl;
//                             // }
//                         }
//                         else
//                         {
//                             nLargeDistancePoints++;
//                             if (nLargeDistancePoints <= 5)
//                             {
//                                 WarningInFunction
//                                     << "固相patch点" << solidGlobalId 
//                                     << "到最近气相点的距离过大: " 
//                                     << minDistance << " m" << nl
//                                     << "固相点坐标: " << solidPt << nl
//                                     << "容差: " << tolerableDistance << " m" << endl;
//                             }
//                         }
//                     }
//                     else
//                     {
//                         if (nLargeDistancePoints <= 5)
//                         {
//                             WarningInFunction
//                                 << "固相patch点" << solidGlobalId 
//                                 << "未找到可用的气相点进行映射" << endl;
//                         }
//                         nLargeDistancePoints++;
//                     }
//                 }               
//                 // 5. 输出映射统计信息
//                 Info<< nl << "基于patch点遍历的映射统计:" << nl;
//                 Info<< "  固相patch总点数: " << solidPatchPoints.size() << nl;
//                 Info<< "  气相patch总点数: " << gasPatchPoints.size() << nl;
//                 Info<< "  处理的固相点数: " << nMappedPoints << nl;
//                 Info<< "  成功建立映射的点数: " << nValidMappedPoints << nl;
//                 Info<< "  距离过大或无可用点的点数: " << nLargeDistancePoints << nl;               
//                 // 统计剩余未使用的气相点
//                 label nUnusedGasPoints = 0;
//                 forAll(gasPatchPointUsed, i)
//                 {
//                     if (!gasPatchPointUsed[i])
//                     {
//                         nUnusedGasPoints++;
//                     }
//                 }
//                 Info<< "  剩余未使用的气相点数: " << nUnusedGasPoints << nl;              
//                 if (nMappedPoints > 0)
//                 {
//                     Info<< "  平均映射距离: " << totalDistance/nMappedPoints << " m" << nl;
//                 }
//                 Info<< "  最大映射距离: " << maxDistance << " m" << nl;
//                 Info<< "  容差: " << tolerableDistance << " m" << nl;                
//                 // 检查映射质量
//                 const scalar mappingRatio = scalar(nValidMappedPoints) / max(solidPatchPoints.size(), 1);                
//                 if (nValidMappedPoints == 0)
//                 {
//                     WarningInFunction
//                         << "未能建立任何点映射关系!" << nl
//                         << "这可能表明两个patch的几何形状不匹配或距离过大" << nl
//                         << "将跳过点映射，使用原始点索引" << endl;
//                 }
//                 else
//                 {
//                     // 映射质量评估
//                     if (maxDistance < 1e-12)
//                     {
//                         Info<< "  映射质量: 优秀 (完全重合)" << nl;
//                     }
//                     else if (maxDistance < 1e-9)
//                     {
//                         Info<< "  映射质量: 很好 (纳米级精度)" << nl;
//                     }
//                     else if (maxDistance < 1e-6)
//                     {
//                         Info<< "  映射质量: 良好 (微米级精度)" << nl;
//                     }
//                     else if (maxDistance < 1e-3)
//                     {
//                         Info<< "  映射质量: 可接受 (毫米级精度)" << nl;
//                     }
//                     else
//                     {
//                         Info<< "  映射质量: 较差 (距离较大)" << nl;
//                     }                    
//                     // 检查映射的完整性
//                     if (mappingRatio > 0.9)
//                     {
//                         Info<< "  映射完整性: 优秀 (>90%)" << nl;
//                     }
//                     else if (mappingRatio > 0.7)
//                     {
//                         Info<< "  映射完整性: 良好 (>70%)" << nl;
//                     }
//                     else if (mappingRatio > 0.5)
//                     {
//                         Info<< "  映射完整性: 一般 (>50%)" << nl;
//                     }
//                     else
//                     {
//                         Info<< "  映射完整性: 较差 (<50%)" << nl;
//                     }
//                 }                
//                 Info<< "最终建立了 " << outputData.solidToGasPointMapping.size() 
//                     << " 个有效的耦合界面点映射关系" << nl;
//                 Info<< "映射成功率: " << (mappingRatio * 100) << "%" << nl;
//                 Info<< "=========================" << nl << endl;
//             }
//             catch (const std::exception& e)
//             {
//                 WarningInFunction
//                     << "无法获取邻域patch信息: " << e.what() << nl
//                     << "将跳过点映射建立" << endl;
//             }
//             catch (...)
//             {
//                 WarningInFunction
//                     << "发生未知异常，将跳过点映射建立" << endl;
//             }
//         }
//         else
//         {
//             Info<< nl << "=== 跳过耦合面点映射 ===" << nl;
//             Info<< "固相patch " << solidToGasPatch.name() 
//                 << " 不是映射边界，边界类型: " << fvp.type() << nl
//                 << "将使用原始点索引，不进行点映射转换" << nl
//                 << "===========================" << nl << endl;
//         }
//     }
//     else
//     {
//         WarningInFunction
//             << "无效的solid_to_gas边界ID: " << couplingPatchID_ << nl
//             << "边界总数: " << pbm.size() << endl;
//     }
// }

void Foam::fvMeshTopoChangers::MeltingDataOutput::applyCouplingPointMapping()
{
    if (debugLevel_ >= 1)
    {
        Info<< nl << "=== 应用耦合面点映射 ===" << nl;
        Info<< "耦合点映射数量: " << solidToGasPointMapping.size() << endl;
    }
    
    // 构建耦合点集合用于快速查找
    labelHashSet couplingPointsSet;
    forAllConstIter(Map<label>, solidToGasPointMapping, iter)
    {
        couplingPointsSet.insert(iter.key());
    }
    
    if (debugLevel_ >= 2)
    {
        Info<< "耦合界面点集合大小: " << couplingPointsSet.size() << endl;
        
        // 调试：显示耦合点ID范围
        if (couplingPointsSet.size() > 0)
        {
            labelList couplingPointsList = couplingPointsSet.sortedToc();
            Info<< "耦合点ID范围: " << couplingPointsList.first() 
                << " 到 " << couplingPointsList.last() << endl;
        }
        
        // 调试：显示收集到的面中的点ID范围
        labelHashSet allFacePointsSet;
        
        forAll(ordinaryBoundaryFaces, i)
        {
            forAll(ordinaryBoundaryFaces[i].pointIds, j)
            {
                allFacePointsSet.insert(ordinaryBoundaryFaces[i].pointIds[j]);
            }
        }
        
        forAll(newSolidToGasBoundaryFaces, i)
        {
            forAll(newSolidToGasBoundaryFaces[i].pointIds, j)
            {
                allFacePointsSet.insert(newSolidToGasBoundaryFaces[i].pointIds[j]);
            }
        }
        
        forAll(internalFaces, i)
        {
            forAll(internalFaces[i].pointIds, j)
            {
                allFacePointsSet.insert(internalFaces[i].pointIds[j]);
            }
        }
        
        if (allFacePointsSet.size() > 0)
        {
            labelList allFacePointsList = allFacePointsSet.sortedToc();
            Info<< "面中点ID范围: " << allFacePointsList.first() 
                << " 到 " << allFacePointsList.last() << endl;
            Info<< "面中点总数: " << allFacePointsSet.size() << endl;
        }
        else
        {
            Info<< "警告：收集到的面中没有任何点!" << endl;
        }
    }
    
    // 统计变量
    label nReplacedPoints = 0;
    label nCheckedPoints = 0;
    
    // 1. 替换普通边界面中的点ID
    label nOrdinaryFaceCouplingPoints = 0;
    forAll(ordinaryBoundaryFaces, i)
    {
        BoundaryFaceInfo& face = ordinaryBoundaryFaces[i];
        bool faceHasCouplingPoints = false;
        
        if (debugLevel_ >= 4)
        {
            Info<< nl << "=== MELTER 替换普通边界面 " << i << " 点ID ===" << nl;
            Info<< "  面ID: " << face.faceId << nl;
            Info<< "  边界: " << face.boundaryName << nl;
            Info<< "  owner单元: " << face.ownCellId << nl;
            Info<< "  替换前点数: " << face.pointIds.size() << nl;
            Info<< "  替换前点列表: " << face.pointIds << nl;
        }
        
        forAll(face.pointIds, j)
        {
            nCheckedPoints++;
            label originalPointId = face.pointIds[j];
            
            if (debugLevel_ >= 4)
            {
                Info<< "  检查点" << j << ": 原ID=" << originalPointId;
            }
            
            if (couplingPointsSet.found(originalPointId))
            {
                faceHasCouplingPoints = true;
                if (debugLevel_ >= 4)
                {
                    Info<< " (是耦合点)";
                }
            }
            
            if (solidToGasPointMapping.found(originalPointId))
            {
                label newPointId = solidToGasPointMapping[originalPointId];
                face.pointIds[j] = newPointId;
                nReplacedPoints++;
                
                if (debugLevel_ >= 4)
                {
                    Info<< " -> 新ID=" << newPointId << " (已替换)" << nl;
                }
                else if (debugLevel_ >= 2)
                {
                    Info<< "  普通边界面" << i << "点" << j << ": " 
                        << originalPointId << " -> " << newPointId << endl;
                }
            }
            else if (debugLevel_ >= 4)
            {
                Info<< " (无映射，保持原ID)" << nl;
            }
        }
        
        if (debugLevel_ >= 4)
        {
            Info<< "  替换后点数: " << face.pointIds.size() << nl;
            Info<< "  替换后点列表: " << face.pointIds << nl;
            Info<< "=======================================" << nl;
        }
        
        if (faceHasCouplingPoints)
        {
            nOrdinaryFaceCouplingPoints++;
            if (debugLevel_ >= 2)
            {
                Info<< "  普通边界面" << i << " (ID=" << face.faceId 
                    << ", 边界=" << face.boundaryName << ") 包含耦合点" << endl;
            }
        }
    }
    
    // 2. 替换新solidToGas边界面中的点ID
    label nNewSolidToGasFaceCouplingPoints = 0;
    forAll(newSolidToGasBoundaryFaces, i)
    {
        NewSolidToGasBoundaryFaceInfo& face = newSolidToGasBoundaryFaces[i];
        bool faceHasCouplingPoints = false;
        
        if (debugLevel_ >= 4)
        {
            Info<< nl << "=== MELTER 替换新solidToGas边界面 " << i << " 点ID ===" << nl;
            Info<< "  面ID: " << face.faceId << nl;
            Info<< "  边界: " << face.boundaryName << nl;
            Info<< "  owner单元: " << face.ownCellId << nl;
            Info<< "  替换前点数: " << face.pointIds.size() << nl;
            Info<< "  替换前点列表: " << face.pointIds << nl;
        }
        
        forAll(face.pointIds, j)
        {
            nCheckedPoints++;
            label originalPointId = face.pointIds[j];
            
            if (debugLevel_ >= 4)
            {
                Info<< "  检查点" << j << ": 原ID=" << originalPointId;
            }
            
            if (couplingPointsSet.found(originalPointId))
            {
                faceHasCouplingPoints = true;
                if (debugLevel_ >= 4)
                {
                    Info<< " (是耦合点)";
                }
            }
            
            if (solidToGasPointMapping.found(originalPointId))
            {
                label newPointId = solidToGasPointMapping[originalPointId];
                face.pointIds[j] = newPointId;
                nReplacedPoints++;
                
                if (debugLevel_ >= 4)
                {
                    Info<< " -> 新ID=" << newPointId << " (已替换)" << nl;
                }
                else if (debugLevel_ >= 2)
                {
                    Info<< "  新solidToGas边界面" << i << "点" << j << ": " 
                        << originalPointId << " -> " << newPointId << endl;
                }
            }
            else if (debugLevel_ >= 4)
            {
                Info<< " (无映射，保持原ID)" << nl;
            }
        }
        
        if (debugLevel_ >= 4)
        {
            Info<< "  替换后点数: " << face.pointIds.size() << nl;
            Info<< "  替换后点列表: " << face.pointIds << nl;
            Info<< "=======================================" << nl;
        }
        
        if (faceHasCouplingPoints)
        {
            nNewSolidToGasFaceCouplingPoints++;
            if (debugLevel_ >= 2)
            {
                Info<< "  新solidToGas边界面" << i << " (ID=" << face.faceId 
                    << ") 包含耦合点" << endl;
            }
        }
    }
    
    // 3. 替换内部面中的点ID
    label nInternalFaceCouplingPoints = 0;
    forAll(internalFaces, i)
    {
        InternalFaceInfo& face = internalFaces[i];
        bool faceHasCouplingPoints = false;
        
        if (debugLevel_ >= 4)
        {
            Info<< nl << "=== MELTER 替换内部面 " << i << " 点ID ===" << nl;
            Info<< "  面ID: " << face.faceId << nl;
            Info<< "  owner单元: " << face.ownCellId << nl;
            Info<< "  neighbor单元: " << face.neiCellId << nl;
            Info<< "  替换前点数: " << face.pointIds.size() << nl;
            Info<< "  替换前点列表: " << face.pointIds << nl;
        }
        
        forAll(face.pointIds, j)
        {
            nCheckedPoints++;
            label originalPointId = face.pointIds[j];
            
            if (debugLevel_ >= 4)
            {
                Info<< "  检查点" << j << ": 原ID=" << originalPointId;
            }
            
            if (couplingPointsSet.found(originalPointId))
            {
                faceHasCouplingPoints = true;
                if (debugLevel_ >= 4)
                {
                    Info<< " (是耦合点)";
                }
            }
            
            if (solidToGasPointMapping.found(originalPointId))
            {
                label newPointId = solidToGasPointMapping[originalPointId];
                face.pointIds[j] = newPointId;
                nReplacedPoints++;
                
                if (debugLevel_ >= 4)
                {
                    Info<< " -> 新ID=" << newPointId << " (已替换)" << nl;
                }
                else if (debugLevel_ >= 2)
                {
                    Info<< "  内部面" << i << "点" << j << ": " 
                        << originalPointId << " -> " << newPointId << endl;
                }
            }
            else if (debugLevel_ >= 4)
            {
                Info<< " (无映射，保持原ID)" << nl;
            }
        }
        
        if (debugLevel_ >= 4)
        {
            Info<< "  替换后点数: " << face.pointIds.size() << nl;
            Info<< "  替换后点列表: " << face.pointIds << nl;
            Info<< "=======================================" << nl;
        }
        
        if (faceHasCouplingPoints)
        {
            nInternalFaceCouplingPoints++;
            if (debugLevel_ >= 2)
            {
                Info<< "  内部面" << i << " (ID=" << face.faceId 
                    << ", owner=" << face.ownCellId << ", nei=" << face.neiCellId 
                    << ") 包含耦合点" << endl;
            }
        }
    }
    
    if (debugLevel_ >= 1)
    {
        Info<< "耦合面点映射应用完成:" << nl
            << "  检查的总点数: " << nCheckedPoints << nl
            << "  实际替换的点数: " << nReplacedPoints << nl
            << "  包含耦合点的普通边界面数: " << nOrdinaryFaceCouplingPoints << nl
            << "  包含耦合点的新solidToGas边界面数: " << nNewSolidToGasFaceCouplingPoints << nl
            << "  包含耦合点的内部面数: " << nInternalFaceCouplingPoints << endl;
        
        if (nCheckedPoints == 0)
        {
            WarningInFunction
                << "没有检查到任何点！这表明收集到的面数据为空" << nl
                << "  普通边界面数: " << ordinaryBoundaryFaces.size() << nl
                << "  新solidToGas边界面数: " << newSolidToGasBoundaryFaces.size() << nl
                << "  内部面数: " << internalFaces.size() << endl;
        }
        else if (nReplacedPoints == 0 && solidToGasPointMapping.size() > 0)
        {
            Info<< "注意: 虽然建立了耦合点映射，但被移除单元的面中不包含耦合界面点" << nl
                << "这表明被移除的单元不在耦合界面附近，或者点ID范围不匹配" << endl;
        }

        Info<< "=========================" << nl << endl;
    }
}

void Foam::fvMeshTopoChangers::MeltingDataOutput::writeToGlobalFile(const polyMesh& mesh) const
{
    const Time& runTime = mesh.time();
    fileName globalDataDir = runTime.rootPath()/runTime.globalCaseName()/"meltingData";
    mkDir(globalDataDir);
    
    fileName dataFile = globalDataDir/("meltingData_" + runTime.name() + ".dat");
    
    OFstream os(dataFile);
    
    if (!os.good())
    {
        FatalErrorInFunction
            << "无法写入融化数据文件: " << dataFile
            << exit(FatalError);
    }
    
    if (debugLevel_ >= 1)
    {
        Info<< "开始写入融化数据文件: " << dataFile << endl;
    }
    
    // 使用简化格式 - 每行一个数据项
    os << "// OpenFOAM Melting Data File - Original IDs Format" << nl
       << "// Time: " << runTime.name() << nl
       << "// All data uses original indices (no reordering)" << nl
       << "// Coupling point mapping already applied" << nl << nl;
    
    // 写入点数据（使用原始ID）
    os << "points" << nl;
    os << removedCellPoints.size() << nl;
    forAll(removedCellPoints, i)
    {
        const PointInfo& pt = removedCellPoints[i];
        os << pt.id << " " << pt.coord.x() << " " << pt.coord.y() << " " << pt.coord.z() << nl;
    }
    os << nl;
    
    // 写入单元数据（新增中心坐标）
    os << "cells" << nl;
    os << removedCells.size() << nl;
    forAll(removedCells, i)
    {
        const CellInfo& cell = removedCells[i];
        // 格式：cellId centroid.x centroid.y centroid.z
        os << cell.cellId << " " 
           << cell.centroid.x() << " " << cell.centroid.y() << " " << cell.centroid.z() << nl;
    }
    os << nl;
    
    // 写入solidToGas边界面信息
    os << "solidToGasBoundaryFaces" << nl;
    os << solidToGasBoundaryFaceIds.size() << nl;
    forAll(solidToGasBoundaryFaceIds, i)
    {
        os << solidToGasBoundaryFaceIds[i] << " "
           << (solidToGasFaceOwnerRemoved[i] ? 1 : 0) << " "
           << solidToGasFaceRemovedCellIds[i] << nl;
    }
    os << nl;
    
    // 写入普通边界面（点ID已经过耦合映射处理）
    os << "ordinaryBoundaryFaces" << nl;
    os << ordinaryBoundaryFaces.size() << nl;
    forAll(ordinaryBoundaryFaces, i)
    {
        const BoundaryFaceInfo& face = ordinaryBoundaryFaces[i];
        os << face.faceId << " " << face.ownCellId << " " << face.boundaryName << " " << face.faceType;
        forAll(face.pointIds, j)
        {
            os << " " << face.pointIds[j];
        }
        os << nl;
    }
    os << nl;
    
    // 写入新solidToGas边界面（点ID已经过耦合映射处理）
    os << "newSolidToGasBoundaryFaces" << nl;
    os << newSolidToGasBoundaryFaces.size() << nl;
    forAll(newSolidToGasBoundaryFaces, i)
    {
        const NewSolidToGasBoundaryFaceInfo& face = newSolidToGasBoundaryFaces[i];
        os << face.faceId << " " << face.ownCellId << " " << face.boundaryName << " " << face.faceType;
        forAll(face.pointIds, j)
        {
            os << " " << face.pointIds[j];
        }
        os << nl;
    }
    os << nl;
    
    // 写入内部面（点ID已经过耦合映射处理）
    os << "internalFaces" << nl;
    os << internalFaces.size() << nl;
    forAll(internalFaces, i)
    {
        const InternalFaceInfo& face = internalFaces[i];
        os << face.faceId << " " << face.ownCellId << " " << face.neiCellId << " " << face.faceType;
        forAll(face.pointIds, j)
        {
            os << " " << face.pointIds[j];
        }
        os << nl;
    }
    
    if (debugLevel_ >= 1)
    {
        Info<< "融化数据已写入全局文件: " << dataFile << nl
            << "包含数据:" << nl
            << "  点数: " << removedCellPoints.size() << nl
            << "  单元数: " << removedCells.size() << nl
            << "  solidToGas边界面数: " << solidToGasBoundaryFaceIds.size() << nl
            << "  普通边界面数: " << ordinaryBoundaryFaces.size() << nl
            << "  新solidToGas边界面数: " << newSolidToGasBoundaryFaces.size() << nl
            << "  内部面数: " << internalFaces.size() << endl;
    }
}

void Foam::fvMeshTopoChangers::melter::debugOutput(const MeltingDataOutput& outputData, label level) const
{
    if (level < 1) return;
    
    Info<< nl << "=== MELTER调试输出 (级别" << level << ") ===" << nl;
    Info<< "时间步: " << mesh().time().timeIndex() << " 时间: " << mesh().time().name() << nl;
    
    // 基本统计信息
    Info<< "数据统计:" << nl
        << "  移除单元数: " << outputData.removedCells.size() << nl
        << "  移除单元的点数: " << outputData.removedCellPoints.size() << nl
        << "  solidToGas边界面数: " << outputData.solidToGasBoundaryFaceIds.size() << nl
        << "  普通边界面数: " << outputData.ordinaryBoundaryFaces.size() << nl
        << "  新solidToGas边界面数: " << outputData.newSolidToGasBoundaryFaces.size() << nl
        << "  内部面数: " << outputData.internalFaces.size() << nl
        << "  耦合点映射数: " << outputData.solidToGasPointMapping.size() << endl;
    
    if (level >= 2)
    {
        // 详细的单元信息
        Info<< nl << "移除的单元详情:" << nl;
        forAll(outputData.removedCells, i)
        {
            Info<< "  单元" << i << ": ID=" << outputData.removedCells[i].cellId << endl;
            if (i >= 10 && outputData.removedCells.size() > 15)
            {
                Info<< "  ... (省略" << (outputData.removedCells.size()-10) << "个单元)" << endl;
                break;
            }
        }
        
        // 点信息
        Info<< nl << "移除单元的点详情 (前10个):" << nl;
        forAll(outputData.removedCellPoints, i)
        {
            const auto& pt = outputData.removedCellPoints[i];
            Info<< "  点" << i << ": ID=" << pt.id << " 坐标=" << pt.coord << endl;
            if (i >= 9) break;
        }
        
        // 耦合面点映射详情
        if (outputData.solidToGasPointMapping.size() > 0)
        {
            Info<< nl << "耦合面点映射详情 (前10个):" << nl;
            label count = 0;
            forAllConstIter(Map<label>, outputData.solidToGasPointMapping, iter)
            {
                Info<< "  固相点" << iter.key() << " -> 气相点" << iter() << endl;
                if (++count >= 10) break;
            }
        }
    }
    
    if (level >= 3)
    {
        // 面的详细信息  
        Info<< nl << "普通边界面详情 (前5个):" << nl;
        forAll(outputData.ordinaryBoundaryFaces, i)
        {
            const auto& face = outputData.ordinaryBoundaryFaces[i];
            Info<< "  面" << i << ": ID=" << face.faceId << " owner=" << face.ownCellId 
                << " 边界=" << face.boundaryName << " 点数=" << face.pointIds.size() << endl;
            if (i >= 4) break;
        }
        
        Info<< nl << "新solidToGas边界面详情 (前5个):" << nl;
        forAll(outputData.newSolidToGasBoundaryFaces, i)
        {
            const auto& face = outputData.newSolidToGasBoundaryFaces[i];
            Info<< "  面" << i << ": ID=" << face.faceId << " owner=" << face.ownCellId 
                << " 点数=" << face.pointIds.size() << endl;
            if (i >= 4) break;
        }
        
        Info<< nl << "内部面详情 (前5个):" << nl;
        forAll(outputData.internalFaces, i)
        {
            const auto& face = outputData.internalFaces[i];
            Info<< "  面" << i << ": ID=" << face.faceId << " owner=" << face.ownCellId 
                << " neighbor=" << face.neiCellId << " 点数=" << face.pointIds.size() << endl;
            if (i >= 4) break;
        }
    }
    
    if (debugLevel_ >= 4)
    {
        Info<< nl << "现有网格统计信息:" << nl;
        Info<< "  网格点数: " << mesh().nPoints() << nl;
        Info<< "  网格单元数: " << mesh().nCells() << nl;
        Info<< "  网格面数: " << mesh().nFaces() << nl;
        Info<< "  内部面数: " << mesh().nInternalFaces() << nl;
        
        const polyBoundaryMesh& pbm = mesh().boundaryMesh();
        Info<< "  边界patch数: " << pbm.size() << nl;
        forAll(pbm, patchi)
        {
            Info<< "    patch" << patchi << ": " << pbm[patchi].name() 
                << " (面数=" << pbm[patchi].size() << ")" << nl;
        }
        
        // 显示耦合patch的详细信息
        if (couplingPatchID_ >= 0 && couplingPatchID_ < pbm.size())
        {
            const polyPatch& couplingPatch = pbm[couplingPatchID_];
            Info<< nl << "耦合patch详细信息:" << nl;
            Info<< "  名称: " << couplingPatch.name() << nl;
            Info<< "  类型: " << couplingPatch.type() << nl;
            Info<< "  面数: " << couplingPatch.size() << nl;
            Info<< "  起始面ID: " << couplingPatch.start() << nl;
            
            // 显示前几个耦合面的信息
            const pointField& meshPoints = mesh().points();
            const faceList& meshFaces = mesh().faces();
            const labelList& faceOwner = mesh().faceOwner();
            
            // 显示前5个耦合面详情
            const label nFacesToShow = min(5, couplingPatch.size());
            Info<< "  前" << nFacesToShow << "个耦合面详情:" << nl;
            for (label i = 0; i < nFacesToShow; i++)
            {
                const label globalFaceI = couplingPatch.start() + i;
                const face& f = meshFaces[globalFaceI];
                const label owner = faceOwner[globalFaceI];
                
                Info<< "    面" << i << " (全局ID=" << globalFaceI << "):" << nl;
                Info<< "      owner单元: " << owner << nl;
                Info<< "      点列表: " << f << nl;
                
                point faceCentre = f.centre(meshPoints);
                vector faceArea = f.area(meshPoints);
                
                Info<< "      面中心: " << faceCentre << nl;
                Info<< "      面法向量: " << faceArea << nl;
                Info<< "      面积: " << mag(faceArea) << nl;
            }
            
            // 显示后5个耦合面详情
            if (couplingPatch.size() > 5)
            {
                const label startIdx = max(0, couplingPatch.size() - 5);
                const label endIdx = couplingPatch.size();
                
                Info<< "  后" << (endIdx - startIdx) << "个耦合面详情:" << nl;
                for (label i = startIdx; i < endIdx; i++)
                {
                    const label globalFaceI = couplingPatch.start() + i;
                    const face& f = meshFaces[globalFaceI];
                    const label owner = faceOwner[globalFaceI];
                    
                    Info<< "    面" << i << " (全局ID=" << globalFaceI << "):" << nl;
                    Info<< "      owner单元: " << owner << nl;
                    Info<< "      点列表: " << f << nl;
                    
                    point faceCentre = f.centre(meshPoints);
                    vector faceArea = f.area(meshPoints);
                    
                    Info<< "      面中心: " << faceCentre << nl;
                    Info<< "      面法向量: " << faceArea << nl;
                    Info<< "      面积: " << mag(faceArea) << nl;
                }
            }
            
            // 如果面数在6-10之间，显示中间的面信息
            if (couplingPatch.size() > 10)
            {
                Info<< "  ... 中间" << (couplingPatch.size() - 10) << "个面未显示 ..." << nl;
            }
            else if (couplingPatch.size() >= 6 && couplingPatch.size() <= 10)
            {
                // 显示中间的面
                Info<< "  中间耦合面详情:" << nl;
                for (label i = 5; i < couplingPatch.size() - 5; i++)
                {
                    const label globalFaceI = couplingPatch.start() + i;
                    const face& f = meshFaces[globalFaceI];
                    const label owner = faceOwner[globalFaceI];
                    
                    Info<< "    面" << i << " (全局ID=" << globalFaceI << "):" << nl;
                    Info<< "      owner单元: " << owner << nl;
                    Info<< "      点列表: " << f << nl;
                    
                    point faceCentre = f.centre(meshPoints);
                    vector faceArea = f.area(meshPoints);
                    
                    Info<< "      面中心: " << faceCentre << nl;
                    Info<< "      面法向量: " << faceArea << nl;
                    Info<< "      面积: " << mag(faceArea) << nl;
                }
            }
        }
    }

    Info<< "=========================" << nl << endl;
}

void Foam::fvMeshTopoChangers::melter::readDict()
{
    meltingTemperature_ = dict_.lookup<scalar>("meltingTemperature");
    
    changeInterval_ = dict_.lookupOrDefault<label>("changeInterval", 1);
    
    if (changeInterval_ < 0)
    {
        FatalIOErrorInFunction(dict_)
            << "非法的changeInterval值 " << changeInterval_ << nl
            << "changeInterval应大于或等于1"
            << exit(FatalIOError);
    }
    
    temperatureFieldName_ = dict_.lookupOrDefault<word>("temperatureField", "T");
    
    // 查找solid_to_gas边界
    findSolidToGasPatch();  // 修改
}

void Foam::fvMeshTopoChangers::melter::findSolidToGasPatch()  // 修改函数名
{
    word solidToGasPatchName = dict_.lookupOrDefault<word>("solidToGasPatch", "solid_to_gas");  // 修改
    
    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    
    couplingPatchID_ = -1;  // 修改
    
    forAll(patches, patchi)
    {
        if (patches[patchi].name() == solidToGasPatchName)
        {
            couplingPatchID_ = patchi;  // 修改
            break;
        }
    }

    if (couplingPatchID_ == -1)  // 修改
    {
        FatalErrorInFunction
            << "找不到solid_to_gas边界 '" << solidToGasPatchName << "'\n"  // 修改
            << "可用的边界有: " << patches.names()
            << exit(FatalError);
    }
}

Foam::labelList 
Foam::fvMeshTopoChangers::melter::selectCellsToRemove() const
{
    if (!mesh().foundObject<volScalarField>(temperatureFieldName_))
    {
        FatalErrorInFunction
            << "找不到温度场 " << temperatureFieldName_
            << exit(FatalError);
    }
    
    const volScalarField& T = 
        mesh().lookupObject<volScalarField>(temperatureFieldName_);
    
    DynamicList<label> cellsToRemove;
    
    // 遍历单元，检查温度是否高于融化温度
    forAll(T, celli)
    {
        if (T[celli] > meltingTemperature_)
        {
            cellsToRemove.append(celli);
        }
    }
    Info<< "找到 " << cellsToRemove.size()
        << " 个单元需要移除 (温度高于融化温度 "
        << meltingTemperature_ << " K)" << endl;

    return cellsToRemove;
}

// 修改changeMesh函数中的相关代码
Foam::autoPtr<Foam::polyTopoChangeMap> 
Foam::fvMeshTopoChangers::melter::changeMesh(const labelList& cellsToRemove)
{
    MeltingDataOutput outputData;
    outputData.debugLevel_ = debugLevel_; 
    outputData.clear();
    
    const pointField& meshPoints = mesh().points();
    const cellList& meshCells = mesh().cells();
    const faceList& meshFaces = mesh().faces();
    const labelList& faceOwner = mesh().faceOwner();
    const labelList& faceNeighbour = mesh().faceNeighbour();
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const vectorField& cellCentres = mesh().cellCentres();  // 新增：获取单元中心

    // 创建移除单元的HashSet用于快速查找
    labelHashSet removedCellsSet;
    forAll(cellsToRemove, i)
    {
        removedCellsSet.insert(cellsToRemove[i]);
    }
    
    Info<< "开始收集和分析数据..." << endl;
    
    // 1. 收集旧solid_to_gas边界信息
    labelHashSet oldSolidToGasBoundaryPointsSet;  // 修改
    labelHashSet oldSolidToGasBoundaryFacesSet;   // 修改
    
    if (couplingPatchID_ >= 0 && couplingPatchID_ < pbm.size())  // 修改
    {
        const polyPatch& solidToGasPatch = pbm[couplingPatchID_];  // 修改
        const labelList& solidToGasMeshPoints = solidToGasPatch.meshPoints();  // 修改
        
        // 收集solid_to_gas边界点
        forAll(solidToGasMeshPoints, i)  // 修改
        {
            oldSolidToGasBoundaryPointsSet.insert(solidToGasMeshPoints[i]);
            outputData.solidToGasBoundaryPoints.append(solidToGasMeshPoints[i]);  // 修改
        }
        
        // 遍历solid_to_gas边界的所有面
        forAll(solidToGasPatch, patchFacei)
        {
            const label globalFaceId = solidToGasPatch.start() + patchFacei;
            const label ownCell = faceOwner[globalFaceId];
            
            outputData.solidToGasBoundaryFaceIds.append(patchFacei);  // 修改
            oldSolidToGasBoundaryFacesSet.insert(globalFaceId);
            
            bool isOwnerRemoved = removedCellsSet.found(ownCell);
            outputData.solidToGasFaceOwnerRemoved.append(isOwnerRemoved);  // 修改

            // 记录单元ID（无论是否被删除）
            outputData.solidToGasFaceRemovedCellIds.append(ownCell);  // 修改
        }
    }
    
    // 新增：建立耦合点映射关系
    buildCouplingPointMapping(outputData);

    // 2. 收集移除单元信息
    labelHashSet allRemovedCellPointsSet;
    
    forAll(cellsToRemove, i)
    {
        const label cellId = cellsToRemove[i];
        const cell& c = meshCells[cellId];
        
        const labelList cellPoints = c.labels(meshFaces);
        
        CellInfo cellInfo(cellId, cellCentres[cellId]);
        outputData.removedCells.append(cellInfo);
        
        forAll(cellPoints, j)
        {
            allRemovedCellPointsSet.insert(cellPoints[j]);
        }
    }
    
    // 3. 从被移除单元的点中去除旧solid_to_gas边界上的点
    labelList filteredPointIds = allRemovedCellPointsSet.sortedToc();
    
    forAll(filteredPointIds, i)
    {
        const label pointId = filteredPointIds[i];
        
        if (!oldSolidToGasBoundaryPointsSet.found(pointId))  // 修改
        {
            PointInfo ptInfo(pointId, meshPoints[pointId]);
            outputData.removedCellPoints.append(ptInfo);
        }
    }
    
    mesh().preChange();
    
    // 标记要移除的单元
    PackedBoolList cellsToKeep(mesh().nCells(), true);
    forAll(cellsToRemove, i)
    {
        cellsToKeep[cellsToRemove[i]] = false;
    }
    
    // 创建网格拓扑变化引擎
    polyTopoChange meshMod(mesh());

    // 移除要删除的单元
    forAll(cellsToRemove, i)
    {
        meshMod.removeCell(cellsToRemove[i], -1);
    }
    
    // 标记已使用的点
    PackedBoolList pointsInUse(mesh().nPoints(), false);
    
     // 4. 处理内部面并分类收集面信息
    for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
    {
        const label own = faceOwner[facei];
        const label nei = faceNeighbour[facei];
        
        const bool ownKeep = cellsToKeep[own];
        const bool neiKeep = cellsToKeep[nei];
        
        if (ownKeep && neiKeep)
        {
            // 两个单元都保留 - 内部面保持不变
            const face& f = meshFaces[facei];
            forAll(f, fp)
            {
                pointsInUse[f[fp]] = true;
            }
        }
        else if (!ownKeep && !neiKeep)
        {
            // 两个单元都移除 - 移除面
            // 检查是否需要收集此面信息
            if (removedCellsSet.found(own) || removedCellsSet.found(nei))
            {
                // 这是被移除单元的内部面
                const face& f = meshFaces[facei];
                labelList facePoints(f.size());
                forAll(f, pointi)
                {
                    facePoints[pointi] = f[pointi];
                }
                
                if (debugLevel_ >= 4)
                {
                    Info<< nl << "=== MELTER 收集内部面信息 ===" << nl;
                    Info<< "  面ID: " << facei << nl;
                    Info<< "  owner单元: " << own << nl;
                    Info<< "  neighbor单元: " << nei << nl;
                    Info<< "  原始面: " << f << nl;
                    Info<< "  收集的点数: " << facePoints.size() << nl;
                    Info<< "  收集的点列表: " << facePoints << nl;
                    
                    // 显示每个点的坐标
                    Info<< "  点坐标详情:" << nl;
                    forAll(facePoints, j)
                    {
                        label pointId = facePoints[j];
                        if (pointId >= 0 && pointId < meshPoints.size())
                        {
                            Info<< "    点" << j << " (ID=" << pointId 
                                << "): " << meshPoints[pointId] << nl;
                        }
                        else
                        {
                            Info<< "    点" << j << " (ID=" << pointId 
                                << "): 无效ID!" << nl;
                        }
                    }
                    Info<< "===========================" << nl;
                }
                
                InternalFaceInfo internalFaceInfo(facei, facePoints, own, nei);
                outputData.internalFaces.append(internalFaceInfo);
            }
            
            meshMod.removeFace(facei, -1);
        }
        else
        {
            // 一个单元保留，一个单元移除 - 将面转换为边界面
            
            label keepCelli, removedCelli;
            bool flipFace = false;
            
            if (ownKeep)
            {
                keepCelli = own;
                removedCelli = nei;
                flipFace = false;
            }
            else
            {
                keepCelli = nei;
                removedCelli = own;
                flipFace = true;
            }
            
            // 获取原始面
            face f = meshFaces[facei];
            if (flipFace)
            {
                f = f.reverseFace();
            }
            
            // 收集面的点信息
            labelList facePoints(f.size());
            forAll(f, pointi)
            {
                facePoints[pointi] = f[pointi];
            }
            
            if (debugLevel_ >= 4)
            {
                Info<< nl << "=== MELTER 收集新solidToGas边界面信息 ===" << nl;
                Info<< "  面ID: " << facei << nl;
                Info<< "  将成为solidToGas边界面" << nl;
                Info<< "  保留单元: " << keepCelli << nl;
                Info<< "  移除单元: " << removedCelli << nl;
                Info<< "  是否翻转: " << (flipFace ? "是" : "否") << nl;
                Info<< "  处理后面: " << f << nl;
                Info<< "  收集的点数: " << facePoints.size() << nl;
                Info<< "  收集的点列表: " << facePoints << nl;
                
                // 显示每个点的坐标
                Info<< "  点坐标详情:" << nl;
                forAll(facePoints, j)
                {
                    label pointId = facePoints[j];
                    if (pointId >= 0 && pointId < meshPoints.size())
                    {
                        Info<< "    点" << j << " (ID=" << pointId 
                            << "): " << meshPoints[pointId] << nl;
                    }
                    else
                    {
                        Info<< "    点" << j << " (ID=" << pointId 
                            << "): 无效ID!" << nl;
                    }
                }
                Info<< "=========================================" << nl;
            }
            
            // 这将成为新的solid_to_gas边界面
            NewSolidToGasBoundaryFaceInfo newSolidToGasFaceInfo(facei, facePoints, removedCelli);
            outputData.newSolidToGasBoundaryFaces.append(newSolidToGasFaceInfo);
            
            // 添加为solid_to_gas边界面
            meshMod.modifyFace
            (
                f,                  
                facei,              
                keepCelli,          
                -1,                 
                flipFace,           
                couplingPatchID_
            );
            
            // 标记此面的点为使用
            forAll(f, fp)
            {
                pointsInUse[f[fp]] = true;
            }
        }
    }
    
    // 5. 处理边界面并分类收集面信息
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        
        forAll(pp, patchFacei)
        {
            const label facei = pp.start() + patchFacei;
            const label own = faceOwner[facei];
            
            if (cellsToKeep[own])
            {
                // 所有者单元保留 - 保持边界面不变
                const face& f = meshFaces[facei];
                forAll(f, fp)
                {
                    pointsInUse[f[fp]] = true;
                }
            }
            else
            {
                // 所有者单元被移除 - 收集此边界面信息
                if (removedCellsSet.found(own) && !oldSolidToGasBoundaryFacesSet.found(facei))
                {
                    // 这是被移除单元的普通边界面（非原solid_to_gas边界）
                    const face& f = meshFaces[facei];
                    labelList facePoints(f.size());
                    forAll(f, pointi)
                    {
                        facePoints[pointi] = f[pointi];
                    }
                    
                    if (debugLevel_ >= 4)
                    {
                        Info<< nl << "=== MELTER 收集普通边界面信息 ===" << nl;
                        Info<< "  面ID: " << facei << nl;
                        Info<< "  patch: " << pp.name() << nl;
                        Info<< "  owner单元: " << own << nl;
                        Info<< "  原始面: " << f << nl;
                        Info<< "  收集的点数: " << facePoints.size() << nl;
                        Info<< "  收集的点列表: " << facePoints << nl;
                        
                        // 显示每个点的坐标
                        Info<< "  点坐标详情:" << nl;
                        forAll(facePoints, j)
                        {
                            label pointId = facePoints[j];
                            if (pointId >= 0 && pointId < meshPoints.size())
                            {
                                Info<< "    点" << j << " (ID=" << pointId 
                                    << "): " << meshPoints[pointId] << nl;
                            }
                            else
                            {
                                Info<< "    点" << j << " (ID=" << pointId 
                                    << "): 无效ID!" << nl;
                            }
                        }
                        Info<< "====================================" << nl;
                    }
                    
                    BoundaryFaceInfo boundaryFaceInfo(facei, facePoints, own, pp.name());
                    outputData.ordinaryBoundaryFaces.append(boundaryFaceInfo);
                }
                
                // 移除边界面
                meshMod.removeFace(facei, -1);
            }
        }
    }
    
    // 移除未使用的点
    forAll(pointsInUse, pointi)
    {
        if (!pointsInUse[pointi])
        {
            meshMod.removePoint(pointi, -1);
        }
    }
    // 新增：应用耦合点映射（替换面中的点ID）
    outputData.applyCouplingPointMapping();
    Info <<"网格变化前：" << endl;
    // 调试输出
    debugOutput(outputData, debugLevel_);
    
    // 写入全局文件
    outputData.writeToGlobalFile(mesh());
    
    // 执行网格变更
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh(), false);
    
    Info <<"网格变化后：" << endl;
    // 调试输出
    debugOutput(outputData, debugLevel_);

    // 显示移除统计
    const label nRemovedCells = returnReduce(cellsToRemove.size(), sumOp<label>());
    const label nTotalCells = returnReduce(mesh().nCells(), sumOp<label>());
    
    Info<< "移除了 " << nRemovedCells << " 个高温单元 (T > " 
        << meltingTemperature_ << ")" << nl
        << "新网格共有 " << nTotalCells << " 个单元。" << endl;
    
    // 更新网格数据
    mesh().topoChange(map);

    return map;
}

// 修改构造函数中的输出信息
Foam::fvMeshTopoChangers::melter::melter(fvMesh& mesh, const dictionary& dict)
:
    fvMeshTopoChanger(mesh),
    dict_(dict),
    debugLevel_(dict_.lookupOrDefault<label>("debugLevel", 1)),
    meltingTemperature_(0.0),
    changeInterval_(1),
    couplingPatchID_(-1),  // 修改
    temperatureFieldName_("T"),
    changedSinceWrite_(false),
    timeIndex_(-1),
    removedCells_(0)
{
    readDict();

    Info<< "已创建melter网格拓扑变化器:" << nl
        << "  融化温度阈值: " << meltingTemperature_ << nl
        << "  变更间隔: " << changeInterval_ << nl
        << "  温度场: " << temperatureFieldName_ << nl
        << "  调试级别: " << debugLevel_ << nl
        << "  solid_to_gas边界: " << mesh.boundaryMesh()[couplingPatchID_].name()
        << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::melter::~melter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::melter::update()
{
    // 每个时间步只执行一次
    if (timeIndex_ == mesh().time().timeIndex())
    {
        return false;
    }
    
    // 更新时间步索引
    timeIndex_ = mesh().time().timeIndex();
    
    // 如果changeInterval为0，则禁用网格变化
    if (changeInterval_ == 0)
    {
        return false;
    }
    
    // 不在初始时间步执行（时间索引为0时），因为此时网格尚未初始化移动
    if 
    (
        mesh().time().timeIndex() > 0
     && mesh().time().timeIndex() % changeInterval_ == 0
    )
    {
        // 选择要移除的单元
        removedCells_ = selectCellsToRemove();
        
        const label nCellsToRemove = 
            returnReduce(removedCells_.size(), sumOp<label>());
        
        if (nCellsToRemove > 0)
        {
            // 执行网格变更
            autoPtr<polyTopoChangeMap> map = changeMesh(removedCells_);
            
            changedSinceWrite_ = true;
            return true;
        }
    }
    
    return false;
}


void Foam::fvMeshTopoChangers::melter::topoChange(const polyTopoChangeMap& map)
{
    //NotImplemented;
    // 仅需更新移除单元的列表
    if (removedCells_.size() > 0)
    {
        labelList newRemovedCells(0);
        forAll(removedCells_, i)
        {
            label oldCelli = removedCells_[i];
            label newCelli = map.reverseCellMap()[oldCelli];
            if (newCelli >= 0)
            {
                newRemovedCells.append(newCelli);
            }
        }
        removedCells_ = newRemovedCells;
    }
}


void Foam::fvMeshTopoChangers::melter::mapMesh(const polyMeshMap& map)
{
    NotImplemented;
}


void Foam::fvMeshTopoChangers::melter::distribute(const polyDistributionMap& map)
{
    //NotImplemented;
    // 仅需分发移除单元的列表
    if (removedCells_.size() > 0)
    {
        //labelList newRemovedCells;
        map.distributeCellIndices(removedCells_);
        //removedCells_.transfer(newRemovedCells);
    }
}


bool Foam::fvMeshTopoChangers::melter::write(const bool write) const
{
    if (changedSinceWrite_)
    {
        // 保存移除的单元信息（可选）
        // ...
        
        return true;
    }
    
    return true;
}


// ************************************************************************* //