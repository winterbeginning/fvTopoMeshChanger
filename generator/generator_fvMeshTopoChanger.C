#include "generator_fvMeshTopoChanger.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"
#include "polyDistributionMap.H"
#include "IFstream.H"
#include "fileName.H"
#include "OSspecific.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(generator, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, generator, fvMesh);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
bool Foam::fvMeshTopoChangers::ReceivedMeltingData::readFromGlobalFile(const fileName& dataFile)
{
    if (!isFile(dataFile))
    {
        Info<< "融化数据文件不存在: " << dataFile << endl;
        return false;
    }
    
    IFstream is(dataFile);
    
    if (!is.good())
    {
        WarningInFunction << "无法打开融化数据文件: " << dataFile << endl;
        return false;
    }
    
    Info<< "开始读取融化数据文件: " << dataFile << endl;
    
    try
    {
        string line;

        // 读取points标识
        while (is.good())
        {
            is.getLine(line);
            if (line == "points")
            {
                break;
            }
        }
        
        if (line == "points")
        {
            label nPoints;
            is >> nPoints;
            
            Info<< "读取点数据: " << nPoints << " 个点" << endl;
            
            points.setSize(nPoints);
            for (label i = 0; i < nPoints; i++)
            {
                label pointId;
                scalar x, y, z;
                is >> pointId >> x >> y >> z;
                
                points[i].originalId = pointId;
                points[i].coord = point(x, y, z);
            }
        }
        
        // 读取cells标识
        while (is.good())
        {
            is.getLine(line);
            if (line == "cells")
            {
                break;
            }
        }
        
        if (line == "cells")
        {
            label nCells;
            is >> nCells;
            
            Info<< "读取单元数据: " << nCells << " 个单元" << endl;
            
            cells.setSize(nCells);
            for (label i = 0; i < nCells; i++)
            {
                label cellId;
                scalar cx, cy, cz;
                
                // 新格式：cellId centroid.x centroid.y centroid.z
                is >> cellId >> cx >> cy >> cz;
                
                cells[i].originalId = cellId;
                cells[i].centroid = point(cx, cy, cz);
            }
        }
        
        // 读取solidToGas边界面信息
        while (is.good())
        {
            is.getLine(line);
            if (line == "solidToGasBoundaryFaces")
            {
                break;
            }
        }
        
        if (line == "solidToGasBoundaryFaces")
        {
            label nFaces;
            is >> nFaces;
            
            Info<< "读取solidToGas边界面数据: " << nFaces << " 个面" << endl;
            
            topBoundaryInfo.faceIds.setSize(nFaces);
            topBoundaryInfo.ownerRemoved.setSize(nFaces);
            topBoundaryInfo.removedCellIds.setSize(nFaces);
            
            for (label i = 0; i < nFaces; i++)
            {
                label faceId, ownerRemovedInt, cellId;
                is >> faceId >> ownerRemovedInt >> cellId;
                
                topBoundaryInfo.faceIds[i] = faceId;
                topBoundaryInfo.ownerRemoved[i] = (ownerRemovedInt == 1);
                topBoundaryInfo.removedCellIds[i] = cellId;
            }
        }
        
        // 读取普通边界面 - 增强调试信息
        while (is.good())
        {
            is.getLine(line);
            if (line == "ordinaryBoundaryFaces")
            {
                break;
            }
        }
        
        if (line == "ordinaryBoundaryFaces")
        {
            label nFaces;
            is >> nFaces;
            
            // Info<< "读取普通边界面: " << nFaces << " 个面" << endl;
            
            ordinaryBoundaryFaces.setSize(nFaces);
            
            // 消耗掉数字后的换行符
            is.getLine(line);
            
            for (label i = 0; i < nFaces; i++)
            {
                // 读取整行
                string faceLine;
                is.getLine(faceLine);
                
                // 跳过空行
                while (faceLine.empty() && is.good())
                {
                    is.getLine(faceLine);
                }
                
                if (faceLine.empty()) 
                {
                    WarningInFunction << "普通边界面数据行为空，跳过面" << i << endl;
                    continue;
                }
                
                // Info<< "解析普通边界面" << i << "行: '" << faceLine << "'" << endl;
                
                // 使用IStringStream解析行数据
                IStringStream lineStream(faceLine);
                
                label faceId, ownCellId, faceType;
                word boundaryName;

                if (lineStream >> faceId >> ownCellId >> boundaryName >> faceType)
                {
                    ordinaryBoundaryFaces[i].originalId = faceId;
                    ordinaryBoundaryFaces[i].ownCellId = ownCellId;
                    ordinaryBoundaryFaces[i].boundaryName = boundaryName;
                    ordinaryBoundaryFaces[i].faceType = faceType;
                    
                    // Info<< "  基本信息解析: faceId=" << faceId 
                    //     << " ownCellId=" << ownCellId 
                    //     << " boundaryName=" << boundaryName 
                    //     << " faceType=" << faceType << endl;
                    
                    // 读取点ID - 使用 OpenFOAM 安全读取方法
                    ordinaryBoundaryFaces[i].pointIds.clear();
                    DynamicList<label> tempPointIds;
                    
                    // Info<< "  开始读取点ID列表..." << endl;
                    label pointId;
                    while (true)
                    {
                        if (lineStream.read(pointId))
                        {
                            tempPointIds.append(pointId);
                            // Info<< "    读取到点ID: " << pointId << endl;
                        }
                        else
                        {
                            if (lineStream.eof())
                            {
                                break;  // 正常结束
                            }
                            else
                            {
                                // 跳过无效标记
                                string badToken;
                                lineStream >> badToken;
                                WarningInFunction << "点ID解析错误，跳过无效标记: " << badToken << endl;
                            }
                        }
                    }
                    
                    // 复制到最终列表
                    ordinaryBoundaryFaces[i].pointIds = tempPointIds;
                    
                    // Info<< "  面" << i << ": ID=" << faceId << " owner=" << ownCellId 
                    //     << " 边界=" << boundaryName << " 点数=" << tempPointIds.size() 
                    //     << " 点列表=" << tempPointIds << endl;
                    
                    // 验证点数
                    if (tempPointIds.size() != 4)
                    {
                        WarningInFunction
                            << "面" << i << " (ID=" << faceId << ") 点数异常: " 
                            << tempPointIds.size() << " (期望4个)" << endl;
                    }
                }
                else
                {
                    WarningInFunction
                        << "无法解析普通边界面行: " << faceLine << endl;
                }
            }
        }
        
        // 读取新solidToGas边界面 - 使用与普通边界面相同的读取方式
        while (is.good())
        {
            is.getLine(line);
            if (line == "newSolidToGasBoundaryFaces")
            {
                break;
            }
        }
        
        if (line == "newSolidToGasBoundaryFaces")
        {
            label nFaces;
            is >> nFaces;
            
            // Info<< "读取新solidToGas边界面: " << nFaces << " 个面" << endl;
            
            newTopBoundaryFaces.setSize(nFaces);
            
            // 消耗掉数字后的换行符
            is.getLine(line);
            
            for (label i = 0; i < nFaces; i++)
            {
                string faceLine;
                is.getLine(faceLine);
                
                // 跳过空行
                while (faceLine.empty() && is.good())
                {
                    is.getLine(faceLine);
                }
                
                if (faceLine.empty()) 
                {
                    WarningInFunction << "新solidToGas边界面数据行为空，跳过面" << i << endl;
                    continue;
                }
                
                // Info<< "解析新solidToGas边界面" << i << "行: '" << faceLine << "'" << endl;
                
                IStringStream lineStream(faceLine);
                
                label faceId, ownCellId, faceType;
                word boundaryName;
                
                if (lineStream >> faceId >> ownCellId >> boundaryName >> faceType)
                {
                    newTopBoundaryFaces[i].originalId = faceId;
                    newTopBoundaryFaces[i].ownCellId = ownCellId;
                    newTopBoundaryFaces[i].boundaryName = boundaryName;
                    newTopBoundaryFaces[i].faceType = faceType;
                    
                    // Info<< "  基本信息解析: faceId=" << faceId 
                    //     << " ownCellId=" << ownCellId 
                    //     << " boundaryName=" << boundaryName 
                    //     << " faceType=" << faceType << endl;
                    
                    // 读取点ID - 使用与普通边界面相同的方法
                    newTopBoundaryFaces[i].pointIds.clear();
                    DynamicList<label> tempPointIds;
                    
                    // Info<< "  开始读取点ID列表..." << endl;
                    label pointId;
                    while (true)
                    {
                        if (lineStream.read(pointId))
                        {
                            tempPointIds.append(pointId);
                            // Info<< "    读取到点ID: " << pointId << endl;
                        }
                        else
                        {
                            if (lineStream.eof())
                            {
                                break;  // 正常结束
                            }
                            else
                            {
                                // 跳过无效标记
                                string badToken;
                                lineStream >> badToken;
                                WarningInFunction << "点ID解析错误，跳过无效标记: " << badToken << endl;
                            }
                        }
                    }
                    
                    // 复制到最终列表
                    newTopBoundaryFaces[i].pointIds = tempPointIds;
                    
                    // Info<< "  面" << i << ": ID=" << faceId << " owner=" << ownCellId 
                    //     << " 边界=" << boundaryName << " 点数=" << tempPointIds.size() 
                    //     << " 点列表=" << tempPointIds << endl;
                    
                    // 验证点数
                    if (tempPointIds.size() != 4)
                    {
                        WarningInFunction
                            << "新solidToGas面" << i << " (ID=" << faceId << ") 点数异常: " 
                            << tempPointIds.size() << " (期望4个)" << endl;
                    }
                }
                else
                {
                    WarningInFunction
                        << "无法解析新solidToGas边界面行: " << faceLine << endl;
                }
            }
        }
        
        // 读取内部面 - 使用与普通边界面相同的读取方式
        while (is.good())
        {
            is.getLine(line);
            if (line == "internalFaces")
            {
                break;
            }
        }
        
        if (line == "internalFaces")
        {
            label nFaces;
            is >> nFaces;
            
            // Info<< "读取内部面: " << nFaces << " 个面" << endl;
            
            internalFaces.setSize(nFaces);
            
            // 消耗掉数字后的换行符
            is.getLine(line);
            
            for (label i = 0; i < nFaces; i++)
            {
                string faceLine;
                is.getLine(faceLine);
                
                // 跳过空行
                while (faceLine.empty() && is.good())
                {
                    is.getLine(faceLine);
                }
                
                if (faceLine.empty()) 
                {
                    WarningInFunction << "内部面数据行为空，跳过面" << i << endl;
                    continue;
                }
                
                // Info<< "解析内部面" << i << "行: '" << faceLine << "'" << endl;
                
                IStringStream lineStream(faceLine);
                
                label faceId, ownCellId, neiCellId, faceType;
                
                if (lineStream >> faceId >> ownCellId >> neiCellId >> faceType)
                {
                    internalFaces[i].originalId = faceId;
                    internalFaces[i].ownCellId = ownCellId;
                    internalFaces[i].neiCellId = neiCellId;
                    internalFaces[i].faceType = faceType;
                    
                    // Info<< "  基本信息解析: faceId=" << faceId 
                    //     << " ownCellId=" << ownCellId 
                    //     << " neiCellId=" << neiCellId 
                    //     << " faceType=" << faceType << endl;
                    
                    // 读取点ID - 使用与普通边界面相同的方法
                    internalFaces[i].pointIds.clear();
                    DynamicList<label> tempPointIds;
                    
                    // Info<< "  开始读取点ID列表..." << endl;
                    label pointId;
                    while (true)
                    {
                        if (lineStream.read(pointId))
                        {
                            tempPointIds.append(pointId);
                            // Info<< "    读取到点ID: " << pointId << endl;
                        }
                        else
                        {
                            if (lineStream.eof())
                            {
                                break;  // 正常结束
                            }
                            else
                            {
                                // 跳过无效标记
                                string badToken;
                                lineStream >> badToken;
                                WarningInFunction << "点ID解析错误，跳过无效标记: " << badToken << endl;
                            }
                        }
                    }
                    
                    // 复制到最终列表
                    internalFaces[i].pointIds = tempPointIds;
                    
                    // Info<< "  面" << i << ": ID=" << faceId << " owner=" << ownCellId 
                    //     << " nei=" << neiCellId << " 点数=" << tempPointIds.size() 
                    //     << " 点列表=" << tempPointIds << endl;
                    
                    // 验证点数
                    if (tempPointIds.size() != 4)
                    {
                        WarningInFunction
                            << "内部面" << i << " (ID=" << faceId << ") 点数异常: " 
                            << tempPointIds.size() << " (期望4个)" << endl;
                    }
                }
                else
                {
                    WarningInFunction
                        << "无法解析内部面行: " << faceLine << endl;
                }
            }
        }
        
        isValid = true;
        sourceTime = dataFile.name(); // 存储文件名作为源时间
        
        Info<< "成功读取融化数据:" << nl
            << "  点数: " << points.size() << nl
            << "  单元数: " << cells.size() << nl
            << "  solidToGas边界面数: " << topBoundaryInfo.faceIds.size() << nl
            << "  普通边界面数: " << ordinaryBoundaryFaces.size() << nl
            << "  新solidToGas边界面数: " << newTopBoundaryFaces.size() << nl
            << "  内部面数: " << internalFaces.size() << endl;
        
        return true;
    }
    catch (const std::exception& e)
    {
        WarningInFunction
            << "读取融化数据时发生错误: " << e.what() << nl
            << "文件: " << dataFile << endl;
        
        // 清理部分读取的数据
        clear();
        return false;
    }
    catch (...)
    {
        WarningInFunction
            << "读取融化数据时发生未知错误" << nl
            << "文件: " << dataFile << endl;
        
        clear();
        return false;
    }
}

void Foam::fvMeshTopoChangers::generator::readDict()
{
    solidRegionName_ = dict_.lookupOrDefault<word>("solidRegion", "solid");
    couplingPatchName_ = dict_.lookupOrDefault<word>("couplingPatch", "gas_to_solid");
    
    // 读取边界名称映射
    readBoundaryMapping();

    // 查找耦合边界
    findCouplingPatch();
}

void Foam::fvMeshTopoChangers::generator::readBoundaryMapping()
{
    boundaryNameMapping_.clear();
    
    if (dict_.found("boundaryMapping"))
    {
        const dictionary& mappingDict = dict_.subDict("boundaryMapping");
        
        Info<< "读取边界名称映射:" << nl;
        
        forAllConstIter(dictionary, mappingDict, iter)
        {
            const word& solidBoundaryName = iter().keyword();
            word gasBoundaryName;
            
            if (iter().isStream())
            {
                ITstream& is = iter().stream();
                is >> gasBoundaryName;
                
                boundaryNameMapping_.insert(solidBoundaryName, gasBoundaryName);
                
                Info<< "  固相边界 '" << solidBoundaryName 
                    << "' -> 气相边界 '" << gasBoundaryName << "'" << nl;
            }
            else
            {
                WarningInFunction
                    << "边界映射项 '" << solidBoundaryName 
                    << "' 格式不正确，跳过" << endl;
            }
        }
        
        Info<< "共读取 " << boundaryNameMapping_.size() << " 个边界名称映射" << endl;
    }
    else
    {
        Info<< "未找到边界名称映射配置，将使用相同的边界名称" << endl;
    }
}

Foam::label Foam::fvMeshTopoChangers::generator::selectBestMasterCell(
    const point& newCellCentre
) const
{
    const vectorField& meshCellCentres = mesh().cellCentres();
    
    if (debugLevel_ >= 3)
    {
        Info<< "为新单元选择母单元: 新单元中心=" << newCellCentre << endl;
    }
    
    // 策略1：找到距离最近的现有单元
    label bestMasterCell = -1;
    scalar minDistance = GREAT;
    
    // 遍历所有现有单元，找到最近的
    forAll(meshCellCentres, cellI)
    {
        // 排除新创建的单元（避免循环引用）
        if (originalToNewCellMap_.found(cellI))
        {
            continue;  // 这是新创建的单元，跳过
        }
        
        scalar distance = mag(newCellCentre - meshCellCentres[cellI]);
        
        if (distance < minDistance)
        {
            minDistance = distance;
            bestMasterCell = cellI;
        }
    }
    
    if (debugLevel_ >= 3)
    {
        Info<< "  选择的母单元: " << bestMasterCell 
            << " 距离: " << minDistance 
            << " 母单元中心: " << (bestMasterCell >= 0 ? meshCellCentres[bestMasterCell] : point::zero) << endl;
    }
    
    return bestMasterCell;
}

Foam::word Foam::fvMeshTopoChangers::generator::getGasBoundaryName(const word& solidBoundaryName) const
{
    if (boundaryNameMapping_.found(solidBoundaryName))
    {
        return boundaryNameMapping_[solidBoundaryName];
    }
    else
    {
        // 如果没有找到映射，使用原名称
        return solidBoundaryName;
    }
}

void Foam::fvMeshTopoChangers::generator::findCouplingPatch()
{
    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    
    couplingPatchID_ = -1;
    
    forAll(patches, patchi)
    {
        if (patches[patchi].name() == couplingPatchName_)
        {
            couplingPatchID_ = patchi;
            break;
        }
    }

    if (couplingPatchID_ == -1)
    {
        FatalErrorInFunction
            << "找不到耦合边界 '" << couplingPatchName_ << "'\n"
            << "可用的边界有: " << patches.names()
            << exit(FatalError);
    }
}

bool Foam::fvMeshTopoChangers::generator::readMeltingData()
{
    const Time& runTime = mesh().time();
    
    // 构建数据文件路径
    fileName globalDataDir = runTime.rootPath()/runTime.globalCaseName()/"meltingData";
    fileName dataFile = globalDataDir/("meltingData_" + runTime.name() + ".dat");
    
    meltingData_.clear();
    
    if (meltingData_.readFromGlobalFile(dataFile))
    {
        meltingData_.sourceTime = runTime.name();
        Info<< "成功读取固相融化数据: " << dataFile << endl;
        
        // 调试级别4：写出读取到的数据进行验证
        if (debugLevel_ >= 4)
        {
            Info<< nl << "=== GENERATOR 读取数据验证 (调试级别4) ===" << nl;
            
            // 验证点数据
            Info<< "读取的点数据验证:" << nl;
            Info<< "  点总数: " << meltingData_.points.size() << nl;
            forAll(meltingData_.points, i)
            {
                const ReceivedPointInfo& pt = meltingData_.points[i];
                Info<< "  点" << i << ": 原ID=" << pt.originalId 
                    << " 坐标=" << pt.coord << nl;
            }
            
            // 验证单元数据
            Info<< nl << "读取的单元数据验证:" << nl;
            Info<< "  单元总数: " << meltingData_.cells.size() << nl;
            forAll(meltingData_.cells, i)
            {
                const ReceivedCellInfo& cell = meltingData_.cells[i];
                Info<< "  单元" << i << ": 原ID=" << cell.originalId << nl;
            }
            
            // 验证普通边界面数据
            Info<< nl << "读取的普通边界面数据验证:" << nl;
            Info<< "  普通边界面总数: " << meltingData_.ordinaryBoundaryFaces.size() << nl;
            forAll(meltingData_.ordinaryBoundaryFaces, i)
            {
                const ReceivedFaceInfo& face = meltingData_.ordinaryBoundaryFaces[i];
                Info<< "  面" << i << ": 原ID=" << face.originalId 
                    << " owner=" << face.ownCellId 
                    << " 边界=" << face.boundaryName 
                    << " 类型=" << face.faceType
                    << " 点数=" << face.pointIds.size()
                    << " 点列表=" << face.pointIds << nl;
            }
            
            // 验证新solidToGas边界面数据
            Info<< nl << "读取的新solidToGas边界面数据验证:" << nl;
            Info<< "  新solidToGas边界面总数: " << meltingData_.newTopBoundaryFaces.size() << nl;
            forAll(meltingData_.newTopBoundaryFaces, i)
            {
                const ReceivedFaceInfo& face = meltingData_.newTopBoundaryFaces[i];
                Info<< "  面" << i << ": 原ID=" << face.originalId 
                    << " owner=" << face.ownCellId 
                    << " 边界=" << face.boundaryName 
                    << " 类型=" << face.faceType
                    << " 点数=" << face.pointIds.size()
                    << " 点列表=" << face.pointIds << nl;
            }
            
            // 验证内部面数据
            Info<< nl << "读取的内部面数据验证:" << nl;
            Info<< "  内部面总数: " << meltingData_.internalFaces.size() << nl;
            forAll(meltingData_.internalFaces, i)
            {
                const ReceivedFaceInfo& face = meltingData_.internalFaces[i];
                Info<< "  面" << i << ": 原ID=" << face.originalId 
                    << " owner=" << face.ownCellId 
                    << " nei=" << face.neiCellId 
                    << " 类型=" << face.faceType
                    << " 点数=" << face.pointIds.size()
                    << " 点列表=" << face.pointIds << nl;
            }
            
            // 验证solidToGas边界信息
            Info<< nl << "读取的solidToGas边界信息验证:" << nl;
            Info<< "  边界面ID数: " << meltingData_.topBoundaryInfo.faceIds.size() << nl;
            Info<< "  ownerRemoved数: " << meltingData_.topBoundaryInfo.ownerRemoved.size() << nl;
            Info<< "  removedCellIds数: " << meltingData_.topBoundaryInfo.removedCellIds.size() << nl;
            
            for (label i = 0; i < min(5, meltingData_.topBoundaryInfo.faceIds.size()); i++)
            {
                Info<< "  边界面" << i << ": faceID=" << meltingData_.topBoundaryInfo.faceIds[i]
                    << " ownerRemoved=" << (meltingData_.topBoundaryInfo.ownerRemoved[i] ? "是" : "否")
                    << " removedCellId=" << meltingData_.topBoundaryInfo.removedCellIds[i] << nl;
            }
            
            Info<< "=== GENERATOR 数据验证完成 ===" << nl << endl;
        }
        
        return true;
    }
    
    return false;
}

void Foam::fvMeshTopoChangers::generator::createPoints(polyTopoChange& meshMod)
{
    if (debugLevel_ >= 1)
    {
        Info<< "创建点: " << meltingData_.points.size() << " 个新点" << endl;
    }
    
    forAll(meltingData_.points, i)
    {
        ReceivedPointInfo& ptInfo = meltingData_.points[i];
        
        // 添加点，接受add函数返回的实际ID
        label newPointId = meshMod.addPoint
        (
            ptInfo.coord,
            -1,
            true
        );
        
        // 建立映射关系：原始ID -> 新ID
        originalToNewPointMap_.insert(ptInfo.originalId, newPointId);
        ptInfo.newId = newPointId;
        
        if (debugLevel_ >= 2)
        {
            Info<< "  点 " << ptInfo.originalId << " -> " << newPointId 
                << " 坐标: " << ptInfo.coord << endl;
        }
    }
    
    if (debugLevel_ >= 1)
    {
        Info<< "点创建完成，建立了 " << originalToNewPointMap_.size() << " 个点映射关系" << endl;
    }
}

void Foam::fvMeshTopoChangers::generator::createCells(polyTopoChange& meshMod)
{
    if (debugLevel_ >= 1)
    {
        Info<< "创建单元: " << meltingData_.cells.size() << " 个新单元" << endl;
    }
    
    forAll(meltingData_.cells, i)
    {
        ReceivedCellInfo& cellInfo = meltingData_.cells[i];
        
        // 关键改进：基于单元中心坐标选择最佳母单元
        label masterCellId = selectBestMasterCell(cellInfo.centroid);

        // 添加单元，接受add函数返回的实际ID
        label newCellId = meshMod.addCell(masterCellId);
        
        // 建立映射关系：原始ID -> 新ID
        originalToNewCellMap_.insert(cellInfo.originalId, newCellId);
        cellInfo.newId = newCellId;
        
        if (debugLevel_ >= 2)
        {
            Info<< "  单元 " << cellInfo.originalId << " -> " << newCellId << endl;
        }
    }
    
    if (debugLevel_ >= 1)
    {
        Info<< "单元创建完成，建立了 " << originalToNewCellMap_.size() << " 个单元映射关系" << endl;
    }
}

// 修复：添加正确的返回类型
Foam::label Foam::fvMeshTopoChangers::generator::getNewPointId(label originalPointId) const
{
    // 首先检查是否在新创建的点映射中
    if (originalToNewPointMap_.found(originalPointId))
    {
        return originalToNewPointMap_[originalPointId];
    }
    
    // 如果不在新创建的点映射中，说明是现有网格中的点，返回原ID
    // 注意：这里避免了耦合面上已经映射过的点被重新映射
    return originalPointId;
}

// 修复：添加正确的返回类型  
Foam::label Foam::fvMeshTopoChangers::generator::getNewCellId(label originalCellId) const
{
    // 首先检查是否在新创建的单元映射中
    if (originalToNewCellMap_.found(originalCellId))
    {
        return originalToNewCellMap_[originalCellId];
    }
    
    // 如果不在新创建的单元映射中，返回原ID
    return originalCellId;
}

void Foam::fvMeshTopoChangers::generator::createOrdinaryBoundaryFaces(polyTopoChange& meshMod)
{
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const pointField& meshPoints = mesh().points();
    
    if (debugLevel_ >= 1)
    {
        Info<< "创建普通边界面: " << meltingData_.ordinaryBoundaryFaces.size() << " 个面" << endl;
    }
    
    forAll(meltingData_.ordinaryBoundaryFaces, i)
    {
        ReceivedFaceInfo& faceInfo = meltingData_.ordinaryBoundaryFaces[i];
        
        const word& solidBoundaryName = faceInfo.boundaryName;
        const word gasBoundaryName = getGasBoundaryName(solidBoundaryName);

        label patchID = -1;
        forAll(pbm, patchi)
        {
            if (pbm[patchi].name() == gasBoundaryName)
            {
                patchID = patchi;
                break;
            }
        }
        
        if (patchID == -1)
        {
            WarningInFunction
                << "找不到气相边界 '" << gasBoundaryName 
                << "' (对应固相边界 '" << solidBoundaryName << "')" << endl;
            continue;
        }
        
        // 调试级别4：输出面创建前的详细信息
        if (debugLevel_ >= 4)
        {
            Info<< nl << "=== 创建普通边界面 " << i << " 详细信息 ===" << nl;
            Info<< "  原始面ID: " << faceInfo.originalId << nl;
            Info<< "  所属边界: " << solidBoundaryName << " -> " << gasBoundaryName << nl;
            Info<< "  目标patch ID: " << patchID << nl;
            Info<< "  原始owner单元ID: " << faceInfo.ownCellId << nl;
            Info<< "  点数: " << faceInfo.pointIds.size() << nl;
            
            Info<< "  原始点ID和坐标:" << nl;
            forAll(faceInfo.pointIds, j)
            {
                label originalPointId = faceInfo.pointIds[j];
                // 查找该点在接收数据中的坐标
                point originalCoord = point::zero;
                bool found = false;
                forAll(meltingData_.points, k)
                {
                    if (meltingData_.points[k].originalId == originalPointId)
                    {
                        originalCoord = meltingData_.points[k].coord;
                        found = true;
                        break;
                    }
                }
                
                Info<< "    点" << j << ": 原ID=" << originalPointId;
                if (found)
                {
                    Info<< " 原坐标=" << originalCoord;
                }
                else
                {
                    Info<< " (坐标未找到)";
                }
                
                label newPointId = getNewPointId(originalPointId);
                Info<< " 新ID=" << newPointId;
                
                if (originalToNewPointMap_.found(originalPointId))
                {
                    // 这是新创建的点，使用原坐标
                    Info<< " 新坐标=" << originalCoord << " (新创建)";
                }
                else
                {
                    // 这是现有网格中的点
                    if (newPointId < meshPoints.size())
                    {
                        Info<< " 新坐标=" << meshPoints[newPointId] << " (现有)";
                    }
                    else
                    {
                        Info<< " (坐标超出范围)";
                    }
                }
                Info<< nl;
            }
        }
        
        // 构建面的点列表，使用新的映射方法
        face f(faceInfo.pointIds.size());
        forAll(faceInfo.pointIds, j)
        {
            f[j] = getNewPointId(faceInfo.pointIds[j]);
        }
        
        // 获取owner单元ID
        label ownCellId = getNewCellId(faceInfo.ownCellId);
        
        // 调试级别4：输出面创建时的最终信息
        if (debugLevel_ >= 4)
        {
            Info<< "  最终面数据:" << nl;
            Info<< "    面点列表: " << f << nl;
            Info<< "    owner单元: " << faceInfo.ownCellId << " -> " << ownCellId << nl;
            
            // 计算面中心和法向量
            point faceCentre = f.centre(meshPoints);
            vector faceNormal = f.area(meshPoints);
            
            Info<< "    面中心: " << faceCentre << nl;
            Info<< "    面法向量: " << faceNormal << nl;
            Info<< "    面积: " << mag(faceNormal) << nl;
        }
        
        const polyPatch& patch = pbm[patchID];
        label referenceFaceID = patch.start();

        // 添加边界面
        label newFaceId = meshMod.addFace
        (
            f,
            ownCellId,
            -1,
            referenceFaceID, // 使用patch的起始面作为参考
            false,
            patchID
        );
        
        faceInfo.newId = newFaceId;
        
        if (debugLevel_ >= 2)
        {
            Info<< "  普通边界面 " << faceInfo.originalId << " -> " << newFaceId
                << " 边界: " << solidBoundaryName << "->" << gasBoundaryName
                << " owner: " << faceInfo.ownCellId << "->" << ownCellId << endl;
        }
        
        if (debugLevel_ >= 4)
        {
            Info<< "  创建成功，新面ID: " << newFaceId << nl;
            Info<< "===========================================" << nl;
        }
    }
}

void Foam::fvMeshTopoChangers::generator::convertTemporaryFacesToCouplingFaces(polyTopoChange& meshMod)
{
    if (debugLevel_ >= 1)
    {
        Info<< "将 " << temporaryBoundaryFaceIds_.size() << " 个临时边界面转换为耦合边界面" << endl;
    }
    
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const faceList& meshFaces = mesh().faces();
    const labelList& faceOwner = mesh().faceOwner();
    const pointField& meshPoints = mesh().points();
    
    // 获取临时patch信息
    label temporaryPatchID = 0;
    const polyPatch& temporaryPatch = pbm[temporaryPatchID];
    
    if (debugLevel_ >= 2)
    {
        Info<< "临时patch信息: " << temporaryPatch.name() 
            << " 起始面ID=" << temporaryPatch.start()
            << " 面数=" << temporaryPatch.size() << endl;
    }

    // 找到所有需要转换的面在临时patch中的位置
    labelList facesToConvert;
    DynamicList<label> tempList;
    
    forAll(temporaryBoundaryFaceIds_, i)
    {
        label faceId = temporaryBoundaryFaceIds_[i];
        
        //检查该面是否在临时patch中
        if (faceId >= 0 && 
            faceId < temporaryPatch.size())
        {
            tempList.append(faceId);
        }
        else
        {
            WarningInFunction
                << "临时面ID " << faceId << " 不在临时patch范围内，跳过" << endl;
        }
    }
    
    facesToConvert = tempList;
    
    if (debugLevel_ >= 2)
    {
        Info<< "找到 " << facesToConvert.size() << " 个有效的临时面需要转换" << endl;
    }
    
    // 按面ID排序，确保从前往后转换
    //Foam::stableSort(facesToConvert);

    // 转换每个临时面为耦合面
    forAll(facesToConvert, i)
    {
        label globalFaceId = temporaryPatch.start() + facesToConvert[facesToConvert.size() - 1 - i]; // 从后往前转换
        const face& f = meshFaces[globalFaceId];
        label owner = faceOwner[globalFaceId];
        
        Info << "正在转换临时面 " << i << ": 全局ID=" << globalFaceId 
            << " owner=" << owner << " 点列表=" << f << endl;

        if (debugLevel_ >= 3)
        {
            Info<< "转换面 " << i << ": 全局ID=" << globalFaceId 
                << " owner=" << owner << " 点列表=" << f << endl;
            
            if (debugLevel_ >= 4)
            {
                point faceCentre = f.centre(meshPoints);
                vector faceArea = f.area(meshPoints);
                Info<< "  面中心: " << faceCentre << " 面积: " << mag(faceArea) << endl;
            }
        }
        // 修改面的patch归属，从临时patch转到耦合patch

        face flippedFace(f);
        flippedFace.flip();
        meshMod.modifyFace
        (
            flippedFace,
            globalFaceId,
            owner,
            -1,        // 边界面没有neighbor
            false,
            couplingPatchID_  // 关键：转换到耦合patch
        );
    }
    
    if (debugLevel_ >= 1)
    {
        Info<< "完成 " << facesToConvert.size() << " 个面的转换" << endl;
    }
}

void Foam::fvMeshTopoChangers::generator::createNewTopBoundaryFacesAsTemporary(polyTopoChange& meshMod)
{
    //const pointField& meshPoints = mesh().points();
    
    if (debugLevel_ >= 1)
    {
        Info<< "创建新top边界面作为临时边界面: " << meltingData_.newTopBoundaryFaces.size() << " 个面" << endl;
    }
    
    // 获取临时边界patch ID（通常是ID=0的patch）
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    label temporaryPatchID = 0; // 使用第一个patch作为临时patch
    
    if (temporaryPatchID >= pbm.size())
    {
        FatalErrorInFunction
            << "临时patch ID " << temporaryPatchID << " 超出范围，共有 " 
            << pbm.size() << " 个patches" << exit(FatalError);
    }
    
    if (debugLevel_ >= 2)
    {
        Info<< "使用临时patch: " << pbm[temporaryPatchID].name() 
            << " (ID=" << temporaryPatchID << ")" << endl;
    }
    
        // 关键修改：先记录临时patch的当前面数，用于计算新面的ID
    const polyPatch& temporaryPatch = pbm[temporaryPatchID];
    label tempPatchStartId = temporaryPatch.start();
    label tempPatchCurrentSize = temporaryPatch.size();
    
    if (debugLevel_ >= 2)
    {
        Info<< "临时patch当前信息: 起始ID=" << tempPatchStartId 
            << " 当前面数=" << tempPatchCurrentSize << endl;
    }

    forAll(meltingData_.newTopBoundaryFaces, i)
    {
        ReceivedFaceInfo& faceInfo = meltingData_.newTopBoundaryFaces[i];
        
        // 调试级别4：输出面创建前的详细信息
        if (debugLevel_ >= 4)
        {
            Info<< nl << "=== 创建临时边界面 " << i << " 详细信息 ===" << nl;
            Info<< "  原始面ID: " << faceInfo.originalId << nl;
            Info<< "  临时patch: " << pbm[temporaryPatchID].name() << " (ID=" << temporaryPatchID << ")" << nl;
            Info<< "  原始owner单元ID: " << faceInfo.ownCellId << nl;
            Info<< "  点数: " << faceInfo.pointIds.size() << nl;
        }
        
        // 构建面的点列表
        face f(faceInfo.pointIds.size());
        forAll(faceInfo.pointIds, j)
        {
            f[j] = getNewPointId(faceInfo.pointIds[j]);
        }
        
        // 获取owner单元ID
        label ownCellId = getNewCellId(faceInfo.ownCellId);
        
        // 添加为临时边界面
        label newFaceId = meshMod.addFace
        (
            f,
            ownCellId,
            -1,
            -1,
            false,
            temporaryPatchID  // 关键：使用临时patch ID
        );
        
        faceInfo.newId = newFaceId;
        
         // 关键修改：计算这个面在changeMesh后在临时patch中的预期位置
        label expectedTempFaceId = tempPatchCurrentSize + i;
        
        // 存储预期的临时边界面ID，而不是实际返回的ID
        temporaryBoundaryFaceIds_.append(expectedTempFaceId);
        
        if (debugLevel_ >= 2)
        {
            Info<< "  临时边界面 " << faceInfo.originalId << " -> " << newFaceId
                << " (临时patch: " << pbm[temporaryPatchID].name() << ")"
                << " owner: " << faceInfo.ownCellId << "->" << ownCellId << endl;
        }
    }
    
    if (debugLevel_ >= 1)
    {
        Info<< "创建了 " << temporaryBoundaryFaceIds_.size() << " 个临时边界面" << endl;
    }
}

void Foam::fvMeshTopoChangers::generator::createInternalFaces(polyTopoChange& meshMod)
{
    if (debugLevel_ >= 1)
    {
        Info<< "创建内部面: " << meltingData_.internalFaces.size() << " 个面" << endl;
    }
    
    forAll(meltingData_.internalFaces, i)
    {
        ReceivedFaceInfo& faceInfo = meltingData_.internalFaces[i];
        
        // 调试级别4：输出面创建前的详细信息
        if (debugLevel_ >= 4)
        {
            Info<< nl << "=== 创建内部面 " << i << " 详细信息 ===" << nl;
            Info<< "  原始面ID: " << faceInfo.originalId << nl;
            Info<< "  原始owner单元ID: " << faceInfo.ownCellId << nl;
            Info<< "  原始neighbor单元ID: " << faceInfo.neiCellId << nl;
            Info<< "  点数: " << faceInfo.pointIds.size() << nl;
            
            Info<< "  原始点ID和坐标:" << nl;
            forAll(faceInfo.pointIds, j)
            {
                label originalPointId = faceInfo.pointIds[j];
                // 查找该点在接收数据中的坐标
                point originalCoord = point::zero;
                bool found = false;
                forAll(meltingData_.points, k)
                {
                    if (meltingData_.points[k].originalId == originalPointId)
                    {
                        originalCoord = meltingData_.points[k].coord;
                        found = true;
                        break;
                    }
                }
                
                Info<< "    点" << j << ": 原ID=" << originalPointId;
                if (found)
                {
                    Info<< " 原坐标=" << originalCoord;
                }
                else
                {
                    Info<< " (坐标未找到)";
                }
                
                label newPointId = getNewPointId(originalPointId);
                Info<< " 新ID=" << newPointId;
                
                if (originalToNewPointMap_.found(originalPointId))
                {
                    // 这是新创建的点，使用原坐标
                    Info<< " 新坐标=" << originalCoord << " (新创建)";
                }
                else
                {
                    // 这是现有网格中的点
                    const pointField& meshPoints = mesh().points();
                    if (newPointId < meshPoints.size())
                    {
                        Info<< " 新坐标=" << meshPoints[newPointId] << " (现有)";
                    }
                    else
                    {
                        Info<< " (坐标超出范围)";
                    }
                }
                Info<< nl;
            }
        }
        
        // 构建面的点列表，使用新的映射方法
        face f(faceInfo.pointIds.size());
        forAll(faceInfo.pointIds, j)
        {
            f[j] = getNewPointId(faceInfo.pointIds[j]);
        }
        
        // 获取owner和neighbor单元ID
        label ownCellId = getNewCellId(faceInfo.ownCellId);
        label neiCellId = getNewCellId(faceInfo.neiCellId);
        
        // 调试级别4：输出面创建时的最终信息
        if (debugLevel_ >= 4)
        {
            Info<< "  最终面数据:" << nl;
            Info<< "    面点列表: " << f << nl;
            Info<< "    owner单元: " << faceInfo.ownCellId << " -> " << ownCellId << nl;
            Info<< "    neighbor单元: " << faceInfo.neiCellId << " -> " << neiCellId << nl;
        }
        
        // 添加内部面
        label newFaceId = meshMod.addFace
        (
            f,
            ownCellId,
            neiCellId,
            -1,
            false,
            -1  // -1表示内部面，无patch
        );
        
        faceInfo.newId = newFaceId;
        
        if (debugLevel_ >= 2)
        {
            Info<< "  内部面 " << faceInfo.originalId << " -> " << newFaceId
                << " owner: " << faceInfo.ownCellId << "->" << ownCellId
                << " nei: " << faceInfo.neiCellId << "->" << neiCellId << endl;
        }
        
        if (debugLevel_ >= 4)
        {
            Info<< "  创建成功，新面ID: " << newFaceId << nl;
            Info<< "===========================================" << nl;
        }
    }
}

void Foam::fvMeshTopoChangers::generator::modifyExistingCouplingFaces(polyTopoChange& meshMod)
{
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const polyPatch& couplingPatch = pbm[couplingPatchID_];
    const labelList& faceOwner = mesh().faceOwner();
    const faceList& meshFaces = mesh().faces();
    const pointField& meshPoints = mesh().points();
    
    if (debugLevel_ >= 1)
    {
        Info<< "修改现有耦合边界面为内部面..." << endl;
    }
    
    forAll(meltingData_.topBoundaryInfo.faceIds, i)
    {
        const label patchFaceI = meltingData_.topBoundaryInfo.faceIds[i];
        const bool ownerRemoved = meltingData_.topBoundaryInfo.ownerRemoved[i];
        const label removedCellId = meltingData_.topBoundaryInfo.removedCellIds[i];
        
        if (patchFaceI < 0 || patchFaceI >= couplingPatch.size())
        {
            WarningInFunction
                << "patch面ID " << patchFaceI << " 超出范围，跳过" << endl;
            continue;
        }
        
        const label globalFaceId = couplingPatch.start() + patchFaceI;
        const label currentOwner = faceOwner[globalFaceId];
        const face& f = meshFaces[globalFaceId];
        
        if (ownerRemoved)
        {
                    // 调试级别4：输出修改前的详细信息
            if (debugLevel_ >= 4)
            {
                Info<< nl << "=== 修改现有耦合面 " << i << " 详细信息 ===" << nl;
                Info<< "  patch面ID: " << patchFaceI << nl;
                Info<< "  全局面ID: " << globalFaceId << nl;
                Info<< "  当前owner单元: " << currentOwner << nl;
                Info<< "  ownerRemoved: " << (ownerRemoved ? "是" : "否") << nl;
                Info<< "  removed单元ID: " << removedCellId << nl;
                
                Info<< "  现有面点信息:" << nl;
                Info<< "    面点列表: " << f << nl;
                forAll(f, j)
                {
                    label pointId = f[j];
                    if (pointId < meshPoints.size())
                    {
                        Info<< "    点" << j << ": ID=" << pointId << " 坐标=" << meshPoints[pointId] << nl;
                    }
                    else
                    {
                        Info<< "    点" << j << ": ID=" << pointId << " (坐标超出范围)" << nl;
                    }
                }
                
                // 计算现有面的几何信息
                point faceCentre = f.centre(meshPoints);
                vector faceArea = f.area(meshPoints);
                
                Info<< "    现有面中心: " << faceCentre << nl;
                Info<< "    现有面法向量: " << faceArea << nl;
                Info<< "    现有面积: " << mag(faceArea) << nl;
            }
            // 获取新的neighbor单元ID
            label newNeiCellId = getNewCellId(removedCellId);
            
            if (debugLevel_ >= 4)
            {
                Info<< "  修改操作:" << nl;
                Info<< "    新neighbor单元: " << removedCellId << " -> " << newNeiCellId << nl;
                Info<< "    修改为内部面: owner=" << currentOwner << " nei=" << newNeiCellId << nl;
            }
            
            if (debugLevel_ >= 2)
            {
                Info<< "  修改耦合面 " << patchFaceI << " (全局ID: " << globalFaceId << ")"
                    << " 为内部面: owner=" << currentOwner << " nei=" 
                    << removedCellId << "->" << newNeiCellId << endl;
            }
            
            meshMod.modifyFace
            (
                f,
                globalFaceId,
                currentOwner,
                newNeiCellId,
                false,
                -1  // patchID=-1表示内部面
            );
            
            if (debugLevel_ >= 4)
            {
                Info<< "  修改成功！" << nl;
                Info<< "===========================================" << nl;
            }
        }
    }
}

void Foam::fvMeshTopoChangers::generator::debugOutput(label level) const
{
    if (level < 1)
    {
        Info<< nl << "现有网格统计信息:" << nl;
        Info<< "  网格点数: " << mesh().nPoints() << nl;
        Info<< "  网格单元数: " << mesh().nCells() << nl;
        Info<< "  网格面数: " << mesh().nFaces() << nl;
        Info<< "  内部面数: " << mesh().nInternalFaces() << nl;
        
        const polyBoundaryMesh& pbm = mesh().boundaryMesh();
        Info<< "  边界patch数: " << pbm.size() << nl;

        return;
    }
    
    Info<< nl << "=== GENERATOR调试输出 (级别" << level << ") ===" << nl;
    Info<< "时间步: " << mesh().time().timeIndex() << " 时间: " << mesh().time().name() << nl;
    
    Info<< "数据统计:" << nl
        << "  接收的点数: " << meltingData_.points.size() << nl
        << "  接收的单元数: " << meltingData_.cells.size() << nl
        << "  普通边界面数: " << meltingData_.ordinaryBoundaryFaces.size() << nl
        << "  新top边界面数: " << meltingData_.newTopBoundaryFaces.size() << nl
        << "  内部面数: " << meltingData_.internalFaces.size() << nl
        << "  需要修改的耦合面数: " << meltingData_.topBoundaryInfo.faceIds.size() << endl;
    
    Info<< "映射关系统计:" << nl
        << "  点映射数: " << originalToNewPointMap_.size() << nl
        << "  单元映射数: " << originalToNewCellMap_.size() << endl;
    
    if (level >= 2)
    {
        Info<< nl << "点映射详情 (前10个):" << nl;
        label count = 0;
        forAllConstIter(Map<label>, originalToNewPointMap_, iter)
        {
            Info<< "  原点" << iter.key() << " -> 新点" << iter() << endl;
            if (++count >= 10) break;
        }
        
        Info<< nl << "单元映射详情 (前10个):" << nl;
        count = 0;
        forAllConstIter(Map<label>, originalToNewCellMap_, iter)
        {
            Info<< "  原单元" << iter.key() << " -> 新单元" << iter() << endl;
            if (++count >= 10) break;
        }
    }
    
    if (level >= 3)
    {
        Info<< nl << "接收的点详情 (前5个):" << nl;
        forAll(meltingData_.points, i)
        {
            const auto& pt = meltingData_.points[i];
            Info<< "  点" << i << ": 原ID=" << pt.originalId << " 新ID=" << pt.newId
                << " 坐标=" << pt.coord << endl;
            if (i >= 4) break;
        }
        
        Info<< nl << "接收的单元详情 (前5个):" << nl;
        forAll(meltingData_.cells, i)
        {
            const auto& cell = meltingData_.cells[i];
            Info<< "  单元" << i << ": 原ID=" << cell.originalId << " 新ID=" << cell.newId << endl;
            if (i >= 4) break;
        }
    }
    
    if (level >= 4)
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

Foam::autoPtr<Foam::polyTopoChangeMap> 
Foam::fvMeshTopoChangers::generator::generateMesh()
{
    // 清空之前的映射
    originalToNewPointMap_.clear();
    originalToNewCellMap_.clear();
    temporaryBoundaryFaceIds_.clear(); // 存储临时边界面ID
    
    mesh().preChange();
    
    polyTopoChange meshMod(mesh());
    
    // 1. 创建点（建立原始ID->新ID映射）
    createPoints(meshMod);
    
    // 2. 创建单元（建立原始ID->新ID映射）
    createCells(meshMod);
    
    // 3. 创建新的top边界面，但设为临时边界面（ID=0的patch）
    createNewTopBoundaryFacesAsTemporary(meshMod);
 
    // 4. 修改现有耦合边界面为内部面
    modifyExistingCouplingFaces(meshMod);

    // 5. 创建普通边界面
    createOrdinaryBoundaryFaces(meshMod);

    // 6. 创建内部面
    createInternalFaces(meshMod);
    
    // 执行第一次网格变更
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh(), false);

    meshChangeStep_ = 1;
      
    mesh().topoChange(map);
    
    mesh().preChange();
    
    polyTopoChange meshMod2(mesh());

    // 将临时边界面转换为耦合边界面
    convertTemporaryFacesToCouplingFaces(meshMod2);

    autoPtr<polyTopoChangeMap> map2 = meshMod2.changeMesh(mesh(), false);

    meshChangeStep_ = 2;

    mesh().topoChange(map2);

    debugOutput(debugLevel_);

    return map2;
}

void Foam::fvMeshTopoChangers::generator::updateCellZones(const polyTopoChangeMap& map)
{
    if (debugLevel_ >= 1)
    {
        Info<< "开始更新cellZones..." << endl;
    }
    
    // 修复1: 从cellMap获取新创建的单元，而不是addedCells()
    const labelList& cellMap = map.cellMap();
    
    // 找到新创建的单元（cellMap中值为-1的表示新创建的单元）
    DynamicList<label> addedCells;
    forAll(cellMap, cellI)
    {
        if (cellMap[cellI] == -1)  // -1表示新创建的单元
        {
            addedCells.append(cellI);
        }
    }
    
    if (addedCells.empty())
    {
        if (debugLevel_ >= 1)
        {
            Info<< "没有新创建的单元，跳过cellZone更新" << endl;
        }
        return;
    }
    
    if (debugLevel_ >= 1)
    {
        Info<< "需要将 " << addedCells.size() << " 个新单元添加到gas cellZone" << endl;
        Info<< "新单元ID列表: " << addedCells << endl;
    }
    
    // 修复2: 获取网格的cellZones，使用正确的类型
    cellZoneList& cellZones = const_cast<cellZoneList&>(mesh().cellZones());
    
    // 查找gas cellZone
    label gasCellZoneID = -1;
    forAll(cellZones, zoneI)
    {
        if (cellZones[zoneI].name() == "gas")
        {
            gasCellZoneID = zoneI;
            break;
        }
    }

    // 获取现有gas cellZone
    cellZone& gasCellZone = cellZones[gasCellZoneID];
    
    if (debugLevel_ >= 2)
    {
        Info<< "找到gas cellZone，当前单元数: " << gasCellZone.size() << endl;
        Info<< "准备添加的新单元: " << addedCells << endl;
    }
    
    // 创建新的单元列表：现有单元 + 新单元
    labelList oldCells = gasCellZone;
    labelList newCells(oldCells.size() + addedCells.size());
    
    // 复制现有单元
    forAll(oldCells, i)
    {
        newCells[i] = oldCells[i];
    }
    
    // 添加新单元
    forAll(addedCells, i)
    {
        newCells[oldCells.size() + i] = addedCells[i];
    }
    
    // 排序（可选，保持有序）
    Foam::stableSort(newCells);
    
    // 修复4: 使用正确的方法更新cellZone
    gasCellZone.transfer(newCells);
    
    if (debugLevel_ >= 1)
    {
        Info<< "已将 " << addedCells.size() << " 个新单元添加到gas cellZone" << endl;
        Info<< "gas cellZone更新后单元数: " << gasCellZone.size() << endl;
    }
    
    if (debugLevel_ >= 2)
    {
        Info<< "更新后的gas cellZone前10个单元: ";
        for (label i = 0; i < min(10, gasCellZone.size()); i++)
        {
            Info<< gasCellZone[i] << " ";
        }
        Info<< endl;
        
        if (gasCellZone.size() > 10)
        {
            Info<< "更新后的gas cellZone后5个单元: ";
            for (label i = max(0, gasCellZone.size()-5); i < gasCellZone.size(); i++)
            {
                Info<< gasCellZone[i] << " ";
            }
            Info<< endl;
        }
    }
    
    // 验证所有cellZones的总单元数
    if (debugLevel_ >= 1)
    {
        label totalCellsInZones = 0;
        forAll(cellZones, zoneI)
        {
            totalCellsInZones += cellZones[zoneI].size();
            Info<< "  cellZone '" << cellZones[zoneI].name() 
                << "': " << cellZones[zoneI].size() << " 个单元" << endl;
        }
        
        Info<< "cellZones总单元数: " << totalCellsInZones 
            << ", 网格总单元数: " << mesh().nCells() << endl;
        
        if (totalCellsInZones != mesh().nCells())
        {
            WarningInFunction
                << "cellZones中的单元数 (" << totalCellsInZones 
                << ") 与网格单元数 (" << mesh().nCells() << ") 不匹配！" << endl;
        }
    }
    
    if (debugLevel_ >= 1)
    {
        Info<< "cellZones更新完成" << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::generator::generator(fvMesh& mesh, const dictionary& dict)
:
    fvMeshTopoChanger(mesh),
    dict_(dict),
    debugLevel_(dict_.lookupOrDefault<label>("debugLevel", 1)), // 新增
    originalToNewPointMap_(),  // 初始化
    originalToNewCellMap_(),   // 初始化
    solidRegionName_("solid"),
    couplingPatchID_(-1),
    couplingPatchName_("gas_to_solid"),
    changedSinceWrite_(false),
    timeIndex_(-1)
{
    readDict();
    
    Info<< "已创建generator网格拓扑变化器:" << nl
        << "  固相区域: " << solidRegionName_ << nl
        << "  耦合边界: " << couplingPatchName_ << nl
        << "  调试级别: " << debugLevel_ << nl
        << "  边界映射数量: " << boundaryNameMapping_.size();
    
    if (boundaryNameMapping_.size() > 0)
    {
        Info<< nl << "  边界映射关系:" << nl;
        
        // 修改：使用正确的迭代器语法
        for 
        (
            HashTable<word, word>::const_iterator iter = boundaryNameMapping_.cbegin();
            iter != boundaryNameMapping_.cend();
            ++iter
        )
        {
            Info<< "    " << iter.key() << " -> " << iter() << nl;
        }
    }
    
    Info<< endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::generator::~generator()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::generator::update()
{
    // 每个时间步只执行一次
    if (timeIndex_ == mesh().time().timeIndex())
    {
        return false;
    }
    
    timeIndex_ = mesh().time().timeIndex();
    
    //尝试读取固相传来的融化数据
    if (readMeltingData() && meltingData_.isValid)
    {
         Info<< "气相区域检测到固相融化数据，开始分两步生成新网格单元..." << endl;
        
        // 第一步：生成网格，新的耦合面创建为临时边界面
        autoPtr<polyTopoChangeMap> map = generateMesh();

        changedSinceWrite_ = true;
        return true;
    }
    
    return false;
}

void Foam::fvMeshTopoChangers::generator::topoChange(const polyTopoChangeMap& map)
{
    if (debugLevel_ >= 1)
    {
        Info<< "generator::topoChange : 处理拓扑变化后的数据清理和cellZone更新" << endl;
    }
    
    switch (meshChangeStep_)
    {
        case 1:  // 第一步完成后
        {
            if (debugLevel_ >= 1)
            {
                Info<< "执行第一步网格变化后的操作..." << endl;
            }
            
            // 第一步后的操作：更新cellZones
            updateCellZones(map);
  
            if (debugLevel_ >= 1)
            {
                Info<< "第一步后操作完成，准备第二步转换" << endl;
            }
            
            break;
        }
        
        case 2:  // 第二步完成后
        {
            if (debugLevel_ >= 1)
            {
                Info<< "执行第二步网格变化后的操作..." << endl;
            }
             
            // 清理所有临时数据
            meltingData_.clear();
            originalToNewPointMap_.clear();
            originalToNewCellMap_.clear();
            temporaryBoundaryFaceIds_.clear();
            
            // 重置步骤标志
            meshChangeStep_ = 0;
            
            if (debugLevel_ >= 1)
            {
                Info<< "两步网格变化全部完成，已清理所有临时数据" << endl;
            }

            break;
        }
        
        default:  // 其他情况（不应该发生）
        {
            WarningInFunction
                << "意外的网格变化步骤: " << meshChangeStep_ << endl;
            break;
        }
    }

    if (debugLevel_ >= 1)
    {
        Info<< "拓扑变化处理完成" << endl;
    }
}

void Foam::fvMeshTopoChangers::generator::mapMesh(const polyMeshMap& map)
{
    NotImplemented;
}

void Foam::fvMeshTopoChangers::generator::distribute(const polyDistributionMap& map)
{
    // 分布式映射处理
    // 如果需要的话，在这里处理并行数据分布
}

bool Foam::fvMeshTopoChangers::generator::write(const bool write) const
{
    return changedSinceWrite_;
}

// ************************************************************************* //