/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
\*---------------------------------------------------------------------------*/

#include "cellRemoval_fvMeshTopoChanger.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"
#include "polyDistributionMap.H" 
#include "OFstream.H"
#include "fileName.H"
#include "OSspecific.H"
#include "mappedFvPatchBaseBase.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(cellRemoval, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, cellRemoval, fvMesh);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshTopoChangers::cellRemoval::RemovalDataOutput::clear()
{
    removedCells.clear();
    removedCellPoints.clear();
    couplingBoundaryPoints.clear();
    couplingBoundaryFaceIds.clear();
    couplingFaceOwnerRemoved.clear();
    couplingFaceRemovedCellIds.clear();
    ordinaryBoundaryFaces.clear();
    newCouplingBoundaryFaces.clear();
    internalFaces.clear();
    couplingPointMapping.clear();
}

void Foam::fvMeshTopoChangers::cellRemoval::RemovalDataOutput::applyCouplingPointMapping()
{
    if (debugLevel_ >= 2)
    {
        Pout<< "应用耦合点映射，映射数量: " << couplingPointMapping.size() << endl;
    }
    
    label nReplacedPoints = 0;
    
    // 替换普通边界面中的点ID
    forAll(ordinaryBoundaryFaces, i)
    {
        BoundaryFaceInfo& face = ordinaryBoundaryFaces[i];
        forAll(face.pointIds, j)
        {
            label originalPointId = face.pointIds[j];
            if (couplingPointMapping.found(originalPointId))
            {
                face.pointIds[j] = couplingPointMapping[originalPointId];
                nReplacedPoints++;
            }
        }
    }
    
    // 替换新耦合边界面中的点ID
    forAll(newCouplingBoundaryFaces, i)
    {
        BoundaryFaceInfo& face = newCouplingBoundaryFaces[i];
        forAll(face.pointIds, j)
        {
            label originalPointId = face.pointIds[j];
            if (couplingPointMapping.found(originalPointId))
            {
                face.pointIds[j] = couplingPointMapping[originalPointId];
                nReplacedPoints++;
            }
        }
    }
    
    // 替换内部面中的点ID
    forAll(internalFaces, i)
    {
        InternalFaceInfo& face = internalFaces[i];
        forAll(face.pointIds, j)
        {
            label originalPointId = face.pointIds[j];
            if (couplingPointMapping.found(originalPointId))
            {
                face.pointIds[j] = couplingPointMapping[originalPointId];
                nReplacedPoints++;
            }
        }
    }
    
    if (debugLevel_ >= 2)
    {
        Pout<< "完成点映射，替换了 " << nReplacedPoints << " 个点" << endl;
    }
}

void Foam::fvMeshTopoChangers::cellRemoval::RemovalDataOutput::writeToFile(const polyMesh& mesh) const
{
    const Time& runTime = mesh.time();
    fileName outputDir;
    
    if (Pstream::parRun())
    {
        // 并行运行：每个处理器在自己的文件夹中写入数据
        outputDir = runTime.rootPath()/runTime.globalCaseName()/"cellRemovalData"/("processor" + Foam::name(Pstream::myProcNo()));
        mkDir(outputDir);
    }
    else
    {
        // 串行运行：按原方式写入
        outputDir = runTime.rootPath()/runTime.globalCaseName()/"cellRemovalData";
        mkDir(outputDir);
    }
    
    fileName dataFile = outputDir/("removalData_" + runTime.name() + ".dat");
    
    OFstream os(dataFile);
    
    if (!os.good())
    {
        FatalErrorInFunction
            << "无法写入单元移除数据文件: " << dataFile
            << exit(FatalError);
    }
    
    // 写入文件头
    os << "// OpenFOAM Cell Removal Data File" << nl
       << "// Time: " << runTime.name() << nl
       << "// Coupling point mapping applied" << nl << nl;
    
    // 写入点数据
    os << "points" << nl << removedCellPoints.size() << nl;
    forAll(removedCellPoints, i)
    {
        const PointInfo& pt = removedCellPoints[i];
        os << pt.id << " " << pt.coord.x() << " " << pt.coord.y() << " " << pt.coord.z() << nl;
    }
    os << nl;
    
    // 写入单元数据
    os << "cells" << nl << removedCells.size() << nl;
    forAll(removedCells, i)
    {
        const CellInfo& cell = removedCells[i];
        os << cell.cellId << " " << cell.centroid.x() << " " << cell.centroid.y() << " " << cell.centroid.z() << nl;
    }
    os << nl;
    
    // 写入耦合边界面信息
    os << "couplingBoundaryFaces" << nl << couplingBoundaryFaceIds.size() << nl;
    forAll(couplingBoundaryFaceIds, i)
    {
        os << couplingBoundaryFaceIds[i] << " "
           << (couplingFaceOwnerRemoved[i] ? 1 : 0) << " "
           << couplingFaceRemovedCellIds[i] << nl;
    }
    os << nl;
    
    // 写入普通边界面
    os << "ordinaryBoundaryFaces" << nl << ordinaryBoundaryFaces.size() << nl;
    forAll(ordinaryBoundaryFaces, i)
    {
        const BoundaryFaceInfo& face = ordinaryBoundaryFaces[i];
        os << face.faceId << " " << face.ownCellId << " " << face.boundaryName;
        forAll(face.pointIds, j)
        {
            os << " " << face.pointIds[j];
        }
        os << nl;
    }
    os << nl;
    
    // 写入新耦合边界面
    os << "newCouplingBoundaryFaces" << nl << newCouplingBoundaryFaces.size() << nl;
    forAll(newCouplingBoundaryFaces, i)
    {
        const BoundaryFaceInfo& face = newCouplingBoundaryFaces[i];
        os << face.faceId << " " << face.ownCellId << " " << face.boundaryName;
        forAll(face.pointIds, j)
        {
            os << " " << face.pointIds[j];
        }
        os << nl;
    }
    os << nl;
    
    // 写入内部面
    os << "internalFaces" << nl << internalFaces.size() << nl;
    forAll(internalFaces, i)
    {
        const InternalFaceInfo& face = internalFaces[i];
        os << face.faceId << " " << face.ownCellId << " " << face.neiCellId;
        forAll(face.pointIds, j)
        {
            os << " " << face.pointIds[j];
        }
        os << nl;
    }
    
    if (debugLevel_ >= 1)
    {
        Pout<< "单元移除数据已写入: " << dataFile << endl;
    }
}

void Foam::fvMeshTopoChangers::cellRemoval::buildCouplingPointMapping(RemovalDataOutput& outputData) const
{
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    
    if (couplingPatchID_ < 0 || couplingPatchID_ >= pbm.size())
    {
        if (debugLevel_ >= 1)
        {
            WarningInFunction << "无效的耦合边界ID: " << couplingPatchID_ << endl;
        }
        return;
    }
    
    const polyPatch& couplingPatch = pbm[couplingPatchID_];
    const fvPatch& fvp = mesh().boundary()[couplingPatchID_];
    
    if (!isA<mappedFvPatchBaseBase>(fvp))
    {
        if (debugLevel_ >= 1)
        {
            Pout<< "耦合边界 " << couplingPatch.name() 
                << " 不是映射边界，跳过点映射建立" << endl;
        }
        return;
    }
    
    const mappedFvPatchBaseBase& mapper = refCast<const mappedFvPatchBaseBase>(fvp);
    
    try
    {
        // 获取邻域patch
        const fvPatch& patchNbr = mapper.nbrFvPatch();
        const polyPatch& polyPatchNbr = patchNbr.patch();
        
        if (debugLevel_ >= 2)
        {
            Pout<< "建立耦合点映射: " << couplingPatch.name() 
                << " <-> " << polyPatchNbr.name() << endl;
        }
        
        // 检查面数量是否匹配
        if (couplingPatch.size() != polyPatchNbr.size())
        {
            WarningInFunction
                << "耦合边界面数量不匹配，跳过点映射建立" << endl;
            return;
        }
        
        if (couplingPatch.size() == 0)
        {
            return;
        }
        
        // 获取网格点坐标
        const pointField& meshPoints1 = mesh().points();
        const pointField& meshPoints2 = patchNbr.boundaryMesh().mesh().points();
        
        // 统计变量
        label nValidMappedPoints = 0;
        scalar totalDistance = 0.0;
        scalar maxDistance = 0.0;
        const scalar tolerableDistance = 1e-3;
        
        // 基于面对应关系建立点映射
        forAll(couplingPatch, facei)
        {
            const face& face1 = couplingPatch[facei];
            const face& face2 = polyPatchNbr[facei];
            
            if (face1.size() != face2.size() || face1.size() == 0)
            {
                continue;
            }
            
            // 提取面点坐标
            pointField facePoints1(face1.size());
            pointField facePoints2(face2.size());
            
            bool validFace = true;
            forAll(face1, pointi)
            {
                label globalId1 = face1[pointi];
                label globalId2 = face2[pointi];
                
                if (globalId1 >= 0 && globalId1 < meshPoints1.size() &&
                    globalId2 >= 0 && globalId2 < meshPoints2.size())
                {
                    facePoints1[pointi] = meshPoints1[globalId1];
                    facePoints2[pointi] = meshPoints2[globalId2];
                }
                else
                {
                    validFace = false;
                    break;
                }
            }
            
            if (!validFace) continue;
            
            // 为每个点找到最近匹配
            forAll(face1, pointi1)
            {
                const point& pt1 = facePoints1[pointi1];
                label globalId1 = face1[pointi1];
                
                label bestMatch = -1;
                scalar minDistance = GREAT;
                
                forAll(face2, pointi2)
                {
                    scalar distance = mag(pt1 - facePoints2[pointi2]);
                    if (distance < minDistance)
                    {
                        minDistance = distance;
                        bestMatch = pointi2;
                    }
                }
                
                if (bestMatch >= 0 && minDistance <= tolerableDistance)
                {
                    label globalId2 = face2[bestMatch];
                    
                    // 检查是否已映射
                    bool alreadyMapped = false;
                    forAllConstIter(Map<label>, outputData.couplingPointMapping, iter)
                    {
                        if (iter() == globalId2)
                        {
                            alreadyMapped = true;
                            break;
                        }
                    }
                    
                    if (!alreadyMapped)
                    {
                        outputData.couplingPointMapping.insert(globalId1, globalId2);
                        nValidMappedPoints++;
                        totalDistance += minDistance;
                        maxDistance = max(maxDistance, minDistance);
                    }
                }
            }
        }
        
        if (debugLevel_ >= 2)
        {
            Pout<< "成功建立 " << nValidMappedPoints << " 个耦合点映射" << endl;
            if (nValidMappedPoints > 0)
            {
                Pout<< "平均距离: " << totalDistance/nValidMappedPoints 
                    << ", 最大距离: " << maxDistance << endl;
            }
        }
    }
    catch (...)
    {
        WarningInFunction << "建立耦合点映射时发生异常" << endl;
    }
}

void Foam::fvMeshTopoChangers::cellRemoval::readDict()
{
    temperatureThreshold_ = dict_.lookup<scalar>("Threshold");
    changeInterval_ = dict_.lookupOrDefault<label>("changeInterval", 1);
    temperatureFieldName_ = dict_.lookupOrDefault<word>("fieldName", "T");
    
    if (changeInterval_ < 0)
    {
        FatalIOErrorInFunction(dict_)
            << "非法的changeInterval值 " << changeInterval_
            << exit(FatalIOError);
    }
    
    findCouplingPatch();
}

void Foam::fvMeshTopoChangers::cellRemoval::findCouplingPatch()
{
    word couplingPatchName = dict_.lookupOrDefault<word>("couplingPatch", "solid_to_gas");
    
    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    
    couplingPatchID_ = -1;
    
    forAll(patches, patchi)
    {
        if (patches[patchi].name() == couplingPatchName)
        {
            couplingPatchID_ = patchi;
            break;
        }
    }

    if (couplingPatchID_ == -1)
    {
        FatalErrorInFunction
            << "找不到耦合边界 '" << couplingPatchName << "'\n"
            << "可用的边界有: " << patches.names()
            << exit(FatalError);
    }
}

Foam::labelList Foam::fvMeshTopoChangers::cellRemoval::selectCellsToRemove() const
{
    if (!mesh().foundObject<volScalarField>(temperatureFieldName_))
    {
        FatalErrorInFunction
            << "找不到温度场 " << temperatureFieldName_
            << exit(FatalError);
    }
    
    const volScalarField& T = mesh().lookupObject<volScalarField>(temperatureFieldName_);
    
    DynamicList<label> cellsToRemove;
    
    forAll(T, celli)
    {
        if (T[celli] > temperatureThreshold_)
        {
            cellsToRemove.append(celli);
        }
    }
    
    if (debugLevel_ >= 1)
    {
        Pout<< "找到 " << cellsToRemove.size() << " 个单元需要移除 (T > " 
            << temperatureThreshold_ << " K)" << endl;
    }

    return cellsToRemove;
}

Foam::autoPtr<Foam::polyTopoChangeMap> 
Foam::fvMeshTopoChangers::cellRemoval::changeMesh(const labelList& cellsToRemove)
{
    RemovalDataOutput outputData;
    outputData.debugLevel_ = debugLevel_;
    outputData.clear();
    
    const pointField& meshPoints = mesh().points();
    const cellList& meshCells = mesh().cells();
    const faceList& meshFaces = mesh().faces();
    const labelList& faceOwner = mesh().faceOwner();
    const labelList& faceNeighbour = mesh().faceNeighbour();
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const vectorField& cellCentres = mesh().cellCentres();

    // 创建移除单元的HashSet
    labelHashSet removedCellsSet;
    forAll(cellsToRemove, i)
    {
        removedCellsSet.insert(cellsToRemove[i]);
    }
    
    // 收集耦合边界信息
    labelHashSet couplingBoundaryPointsSet;
    labelHashSet couplingBoundaryFacesSet;
    
    if (couplingPatchID_ >= 0 && couplingPatchID_ < pbm.size())
    {
        const polyPatch& couplingPatch = pbm[couplingPatchID_];
        const labelList& couplingMeshPoints = couplingPatch.meshPoints();
        
        forAll(couplingMeshPoints, i)
        {
            couplingBoundaryPointsSet.insert(couplingMeshPoints[i]);
            outputData.couplingBoundaryPoints.append(couplingMeshPoints[i]);
        }
        
        forAll(couplingPatch, patchFacei)
        {
            const label globalFaceId = couplingPatch.start() + patchFacei;
            const label ownCell = faceOwner[globalFaceId];
            
            outputData.couplingBoundaryFaceIds.append(patchFacei);
            couplingBoundaryFacesSet.insert(globalFaceId);
            
            bool isOwnerRemoved = removedCellsSet.found(ownCell);
            outputData.couplingFaceOwnerRemoved.append(isOwnerRemoved);
            outputData.couplingFaceRemovedCellIds.append(ownCell);
        }
    }
    
    // 建立耦合点映射
    buildCouplingPointMapping(outputData);

    // 收集移除单元信息
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
    
    // 过滤掉耦合边界上的点
    labelList filteredPointIds = allRemovedCellPointsSet.sortedToc();
    
    forAll(filteredPointIds, i)
    {
        const label pointId = filteredPointIds[i];
        
        if (!couplingBoundaryPointsSet.found(pointId))
        {
            PointInfo ptInfo(pointId, meshPoints[pointId]);
            outputData.removedCellPoints.append(ptInfo);
        }
    }
    
    mesh().preChange();
    
    // 标记要保留的单元
    PackedBoolList cellsToKeep(mesh().nCells(), true);
    forAll(cellsToRemove, i)
    {
        cellsToKeep[cellsToRemove[i]] = false;
    }
    
    // 创建网格拓扑变化引擎
    polyTopoChange meshMod(mesh());

    // 移除单元
    forAll(cellsToRemove, i)
    {
        meshMod.removeCell(cellsToRemove[i], -1);
    }
    
    // 标记使用的点
    PackedBoolList pointsInUse(mesh().nPoints(), false);
    
    // 处理内部面
    for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
    {
        const label own = faceOwner[facei];
        const label nei = faceNeighbour[facei];
        
        const bool ownKeep = cellsToKeep[own];
        const bool neiKeep = cellsToKeep[nei];
        
        if (ownKeep && neiKeep)
        {
            // 保留内部面
            const face& f = meshFaces[facei];
            forAll(f, fp)
            {
                pointsInUse[f[fp]] = true;
            }
        }
        else if (!ownKeep && !neiKeep)
        {
            // 移除内部面，收集信息
            if (removedCellsSet.found(own) || removedCellsSet.found(nei))
            {
                const face& f = meshFaces[facei];
                labelList facePoints(f);
                
                InternalFaceInfo internalFaceInfo(facei, facePoints, own, nei);
                outputData.internalFaces.append(internalFaceInfo);
            }
            
            meshMod.removeFace(facei, -1);
        }
        else
        {
            // 转换为边界面
            label keepCelli = ownKeep ? own : nei;
            label removedCelli = ownKeep ? nei : own;
            bool flipFace = !ownKeep;
            
            face f = meshFaces[facei];
            if (flipFace)
            {
                f = f.reverseFace();
            }
            
            labelList facePoints(f);
            BoundaryFaceInfo newCouplingFaceInfo(facei, facePoints, removedCelli, "coupling");
            outputData.newCouplingBoundaryFaces.append(newCouplingFaceInfo);
            
            meshMod.modifyFace(f, facei, keepCelli, -1, flipFace, couplingPatchID_);
            
            forAll(f, fp)
            {
                pointsInUse[f[fp]] = true;
            }
        }
    }
    
    // 处理边界面
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        
        forAll(pp, patchFacei)
        {
            const label facei = pp.start() + patchFacei;
            const label own = faceOwner[facei];
            
            if (cellsToKeep[own])
            {
                // 保留边界面
                const face& f = meshFaces[facei];
                forAll(f, fp)
                {
                    pointsInUse[f[fp]] = true;
                }
            }
            else
            {
                // 收集边界面信息
                if (removedCellsSet.found(own) && !couplingBoundaryFacesSet.found(facei))
                {
                    const face& f = meshFaces[facei];
                    labelList facePoints(f);
                    
                    BoundaryFaceInfo boundaryFaceInfo(facei, facePoints, own, pp.name());
                    outputData.ordinaryBoundaryFaces.append(boundaryFaceInfo);
                }
                
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
    
    // 应用耦合点映射
    outputData.applyCouplingPointMapping();
    
    // 调试输出和写入文件
    debugOutput(outputData);
    outputData.writeToFile(mesh());
    
    // 执行网格变更
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh(), false);
    
    if (debugLevel_ >= 1)
    {
        Pout<< "移除了 " << cellsToRemove.size() << " 个高温单元，"
        << "新网格共有 " << mesh().nCells() << " 个单元" << endl;
    }

    // 更新网格数据
    mesh().topoChange(map);
    
    return map;
}

void Foam::fvMeshTopoChangers::cellRemoval::debugOutput(const RemovalDataOutput& outputData) const
{
    if (debugLevel_ < 1) return;
    
    Pout<< nl << "=== 单元移除调试信息 ===" << nl;
    Pout<< "时间: " << mesh().time().name() << nl;
    Pout<< "移除单元数: " << outputData.removedCells.size() << nl;
    Pout<< "移除点数: " << outputData.removedCellPoints.size() << nl;
    Pout<< "耦合点映射数: " << outputData.couplingPointMapping.size() << nl;
    Pout<< "普通边界面数: " << outputData.ordinaryBoundaryFaces.size() << nl;
    Pout<< "新耦合边界面数: " << outputData.newCouplingBoundaryFaces.size() << nl;
    Pout<< "内部面数: " << outputData.internalFaces.size() << nl;
    
    if (debugLevel_ >= 2)
    {
        Pout<< nl << "网格统计:" << nl;
        Pout<< "  总点数: " << mesh().nPoints() << nl;
        Pout<< "  总单元数: " << mesh().nCells() << nl;
        Pout<< "  总面数: " << mesh().nFaces() << nl;
        Pout<< "  内部面数: " << mesh().nInternalFaces() << nl;
        
        if (couplingPatchID_ >= 0)
        {
            const polyPatch& couplingPatch = mesh().boundaryMesh()[couplingPatchID_];
            Pout<< "  耦合边界: " << couplingPatch.name() 
                << " (面数=" << couplingPatch.size() << ")" << nl;
        }
    }
    
    Pout<< "=========================" << nl << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::cellRemoval::cellRemoval(fvMesh& mesh, const dictionary& dict)
:
    fvMeshTopoChanger(mesh),
    dict_(dict),
    debugLevel_(dict_.lookupOrDefault<label>("debugLevel", 1)),
    temperatureThreshold_(0.0),
    changeInterval_(1),
    couplingPatchID_(-1),
    temperatureFieldName_("T"),
    changedSinceWrite_(false),
    timeIndex_(-1),
    removedCells_(0)
{
    readDict();

    if (debugLevel_ >= 1)
    {
        Pout<< "创建单元移除网格拓扑变化器:" << nl
            << "  温度阈值: " << temperatureThreshold_ << nl
            << "  变更间隔: " << changeInterval_ << nl
            << "  温度场: " << temperatureFieldName_ << nl
            << "  耦合边界: " << mesh.boundaryMesh()[couplingPatchID_].name()
            << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::cellRemoval::~cellRemoval()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::cellRemoval::update()
{
    // 每个时间步只执行一次
    if (timeIndex_ == mesh().time().timeIndex())
    {
        return false;
    }
    
    timeIndex_ = mesh().time().timeIndex();
    
    if (changeInterval_ == 0)
    {
        return false;
    }
    
    if (mesh().time().timeIndex() > 0 && mesh().time().timeIndex() % changeInterval_ == 0)
    {
        removedCells_ = selectCellsToRemove();
        
        const label nLocalCellsToRemove = removedCells_.size();
        const label nGlobalCellsToRemove = returnReduce(nLocalCellsToRemove, sumOp<label>());

        if (nGlobalCellsToRemove > 0)
        {
            autoPtr<polyTopoChangeMap> map = changeMesh(removedCells_);
            changedSinceWrite_ = true;
            return true;
        }
    }
    
    return false;
}

void Foam::fvMeshTopoChangers::cellRemoval::topoChange(const polyTopoChangeMap& map)
{
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

void Foam::fvMeshTopoChangers::cellRemoval::mapMesh(const polyMeshMap& map)
{
    NotImplemented;
}

void Foam::fvMeshTopoChangers::cellRemoval::distribute(const polyDistributionMap& map)
{
    if (removedCells_.size() > 0)
    {
        map.distributeCellIndices(removedCells_);
    }
}

bool Foam::fvMeshTopoChangers::cellRemoval::write(const bool write) const
{
    return changedSinceWrite_;
}

// ************************************************************************* //