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

#include "cellAddition_fvMeshTopoChanger.H"
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
    defineTypeNameAndDebug(cellAddition, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, cellAddition, fvMesh);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshTopoChangers::cellAddition::AdditionDataInput::clear()
{
    points.clear();
    cells.clear();
    ordinaryBoundaryFaces.clear();
    newCouplingBoundaryFaces.clear();
    internalFaces.clear();
    couplingBoundaryInfo.clear();
    isValid = false;
    sourceTime = word::null;
}

bool Foam::fvMeshTopoChangers::cellAddition::AdditionDataInput::readFromFile(const fileName& dataFile)
{
    if (!isFile(dataFile))
    {
        return false;
    }
    
    IFstream is(dataFile);
    
    if (!is.good())
    {
        WarningInFunction << "无法打开数据文件: " << dataFile << endl;
        return false;
    }
    
    try
    {
        string line;

        // 读取points
        while (is.good())
        {
            is.getLine(line);
            if (line == "points") break;
        }
        
        if (line == "points")
        {
            label nPoints;
            is >> nPoints;
            
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
        
        // 读取cells
        while (is.good())
        {
            is.getLine(line);
            if (line == "cells") break;
        }
        
        if (line == "cells")
        {
            label nCells;
            is >> nCells;
            
            cells.setSize(nCells);
            for (label i = 0; i < nCells; i++)
            {
                label cellId;
                scalar cx, cy, cz;
                is >> cellId >> cx >> cy >> cz;
                
                cells[i].originalId = cellId;
                cells[i].centroid = point(cx, cy, cz);
            }
        }
        
        // 读取耦合边界面信息
        while (is.good())
        {
            is.getLine(line);
            if (line == "couplingBoundaryFaces") break;
        }
        
        if (line == "couplingBoundaryFaces")
        {
            label nFaces;
            is >> nFaces;
            
            couplingBoundaryInfo.faceIds.setSize(nFaces);
            couplingBoundaryInfo.ownerRemoved.setSize(nFaces);
            couplingBoundaryInfo.removedCellIds.setSize(nFaces);
            
            for (label i = 0; i < nFaces; i++)
            {
                label faceId, ownerRemovedInt, cellId;
                is >> faceId >> ownerRemovedInt >> cellId;
                
                couplingBoundaryInfo.faceIds[i] = faceId;
                couplingBoundaryInfo.ownerRemoved[i] = (ownerRemovedInt == 1);
                couplingBoundaryInfo.removedCellIds[i] = cellId;
            }
        }
        
        // 读取普通边界面
        while (is.good())
        {
            is.getLine(line);
            if (line == "ordinaryBoundaryFaces") break;
        }
        
        if (line == "ordinaryBoundaryFaces")
        {
            label nFaces;
            is >> nFaces;
            
            ordinaryBoundaryFaces.setSize(nFaces);
            is.getLine(line); // 消耗换行符
            
            for (label i = 0; i < nFaces; i++)
            {
                string faceLine;
                is.getLine(faceLine);
                
                while (faceLine.empty() && is.good())
                {
                    is.getLine(faceLine);
                }
                
                if (!faceLine.empty())
                {
                    IStringStream lineStream(faceLine);
                    
                    label faceId, ownCellId;
                    word boundaryName;

                    if (lineStream >> faceId >> ownCellId >> boundaryName)
                    {
                        ordinaryBoundaryFaces[i].originalId = faceId;
                        ordinaryBoundaryFaces[i].ownCellId = ownCellId;
                        ordinaryBoundaryFaces[i].boundaryName = boundaryName;
                        
                        // 读取点ID
                        ordinaryBoundaryFaces[i].pointIds.clear();
                        DynamicList<label> tempPointIds;
                        
                        label pointId;
                        while (lineStream.read(pointId))
                        {
                            tempPointIds.append(pointId);
                        }
                        ordinaryBoundaryFaces[i].pointIds = tempPointIds;
                    }
                }
            }
        }
        
        // 读取新耦合边界面
        while (is.good())
        {
            is.getLine(line);
            if (line == "newCouplingBoundaryFaces") break;
        }
        
        if (line == "newCouplingBoundaryFaces")
        {
            label nFaces;
            is >> nFaces;
            
            newCouplingBoundaryFaces.setSize(nFaces);
            is.getLine(line); // 消耗换行符
            
            for (label i = 0; i < nFaces; i++)
            {
                string faceLine;
                is.getLine(faceLine);
                
                while (faceLine.empty() && is.good())
                {
                    is.getLine(faceLine);
                }
                
                if (!faceLine.empty())
                {
                    IStringStream lineStream(faceLine);
                    
                    label faceId, ownCellId;
                    word boundaryName;
                    
                    if (lineStream >> faceId >> ownCellId >> boundaryName)
                    {
                        newCouplingBoundaryFaces[i].originalId = faceId;
                        newCouplingBoundaryFaces[i].ownCellId = ownCellId;
                        newCouplingBoundaryFaces[i].boundaryName = boundaryName;
                        
                        // 读取点ID
                        newCouplingBoundaryFaces[i].pointIds.clear();
                        DynamicList<label> tempPointIds;
                        
                        label pointId;
                        while (lineStream.read(pointId))
                        {
                            tempPointIds.append(pointId);
                        }
                        
                        newCouplingBoundaryFaces[i].pointIds = tempPointIds;
                    }
                }
            }
        }
        
        // 读取内部面
        while (is.good())
        {
            is.getLine(line);
            if (line == "internalFaces") break;
        }
        
        if (line == "internalFaces")
        {
            label nFaces;
            is >> nFaces;
            
            internalFaces.setSize(nFaces);
            is.getLine(line); // 消耗换行符
            
            for (label i = 0; i < nFaces; i++)
            {
                string faceLine;
                is.getLine(faceLine);
                
                while (faceLine.empty() && is.good())
                {
                    is.getLine(faceLine);
                }
                
                if (!faceLine.empty())
                {
                    IStringStream lineStream(faceLine);
                    label faceId, ownCellId, neiCellId;

                    if (lineStream >> faceId >> ownCellId >> neiCellId)
                    {
                        internalFaces[i].originalId = faceId;
                        internalFaces[i].ownCellId = ownCellId;
                        internalFaces[i].neiCellId = neiCellId;
                        
                        // 读取点ID
                        internalFaces[i].pointIds.clear();
                        DynamicList<label> tempPointIds;
                        
                        label pointId;
                        while (lineStream.read(pointId))
                        {
                            tempPointIds.append(pointId);
                        }
                        
                        internalFaces[i].pointIds = tempPointIds;
                    }
                }
            }
        }
        
        isValid = true;
        sourceTime = dataFile.name();
        
        return true;
    }
    catch (...)
    {
        WarningInFunction << "读取数据文件时发生错误: " << dataFile << endl;
        clear();
        return false;
    }
}

void Foam::fvMeshTopoChangers::cellAddition::readDict()
{
    sourceRegionName_ = dict_.lookupOrDefault<word>("sourceRegion", "solid");
    couplingPatchName_ = dict_.lookupOrDefault<word>("couplingPatch", "gas_to_solid");
    
    readBoundaryMapping();
    findCouplingPatch();
}

void Foam::fvMeshTopoChangers::cellAddition::readBoundaryMapping()
{
    boundaryNameMapping_.clear();
    
    if (dict_.found("boundaryMapping"))
    {
        const dictionary& mappingDict = dict_.subDict("boundaryMapping");
        
        if (debugLevel_ >= 1)
        {
            Info<< "读取边界名称映射:" << nl;
        }
        
        forAllConstIter(dictionary, mappingDict, iter)
        {
            const word& sourceBoundaryName = iter().keyword();
            word targetBoundaryName;
            
            if (iter().isStream())
            {
                ITstream& is = iter().stream();
                is >> targetBoundaryName;
                
                boundaryNameMapping_.insert(sourceBoundaryName, targetBoundaryName);
                
                if (debugLevel_ >= 1)
                {
                    Info<< "  " << sourceBoundaryName << " -> " << targetBoundaryName << nl;
                }
            }
        }
        
        if (debugLevel_ >= 1)
        {
            Info<< "共读取 " << boundaryNameMapping_.size() << " 个边界映射" << endl;
        }
    }
}

Foam::word Foam::fvMeshTopoChangers::cellAddition::getTargetBoundaryName(const word& sourceBoundaryName) const
{
    if (boundaryNameMapping_.found(sourceBoundaryName))
    {
        return boundaryNameMapping_[sourceBoundaryName];
    }
    else
    {
        return sourceBoundaryName;
    }
}

void Foam::fvMeshTopoChangers::cellAddition::findCouplingPatch()
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

bool Foam::fvMeshTopoChangers::cellAddition::readAdditionData()
{
    const Time& runTime = mesh().time();
    
    fileName globalDataDir = runTime.rootPath()/runTime.globalCaseName()/"cellRemovalData";
    fileName dataFile = globalDataDir/("removalData_" + runTime.name() + ".dat");
    
    additionData_.clear();
    
    if (additionData_.readFromFile(dataFile))
    {
        additionData_.sourceTime = runTime.name();
        
        if (debugLevel_ >= 1)
        {
            Info<< "成功读取单元添加数据: " << dataFile << nl
                << "  点数: " << additionData_.points.size() << nl
                << "  单元数: " << additionData_.cells.size() << nl
                << "  普通边界面数: " << additionData_.ordinaryBoundaryFaces.size() << nl
                << "  新耦合边界面数: " << additionData_.newCouplingBoundaryFaces.size() << nl
                << "  内部面数: " << additionData_.internalFaces.size() << endl;
        }
        
        return true;
    }
    
    return false;
}

void Foam::fvMeshTopoChangers::cellAddition::createPoints(polyTopoChange& meshMod)
{
    if (debugLevel_ >= 1)
    {
        Info<< "创建点: " << additionData_.points.size() << " 个新点" << endl;
    }
    
    forAll(additionData_.points, i)
    {
        PointInfo& ptInfo = additionData_.points[i];
        
        label newPointId = meshMod.addPoint(ptInfo.coord, -1, true);
        
        originalToNewPointMap_.insert(ptInfo.originalId, newPointId);
        ptInfo.newId = newPointId;
        
        if (debugLevel_ >= 3)
        {
            Info<< "  点 " << ptInfo.originalId << " -> " << newPointId 
                << " 坐标: " << ptInfo.coord << endl;
        }
    }
    
    if (debugLevel_ >= 1)
    {
        Info<< "点创建完成，建立了 " << originalToNewPointMap_.size() << " 个点映射" << endl;
    }
}

void Foam::fvMeshTopoChangers::cellAddition::createCells(polyTopoChange& meshMod)
{
    if (debugLevel_ >= 1)
    {
        Info<< "创建单元: " << additionData_.cells.size() << " 个新单元" << endl;
    }
    
    forAll(additionData_.cells, i)
    {
        CellInfo& cellInfo = additionData_.cells[i];

        label newCellId = meshMod.addCell(-1);
        
        originalToNewCellMap_.insert(cellInfo.originalId, newCellId);
        cellInfo.newId = newCellId;
        
        if (debugLevel_ >= 3)
        {
            Info<< "  单元 " << cellInfo.originalId << " -> " << newCellId << endl;
        }
    }
    
    if (debugLevel_ >= 1)
    {
        Info<< "单元创建完成，建立了 " << originalToNewCellMap_.size() << " 个单元映射" << endl;
    }
}

Foam::label Foam::fvMeshTopoChangers::cellAddition::getNewPointId(label originalPointId) const
{
    if (originalToNewPointMap_.found(originalPointId))
    {
        return originalToNewPointMap_[originalPointId];
    }
    
    return originalPointId;
}

Foam::label Foam::fvMeshTopoChangers::cellAddition::getNewCellId(label originalCellId) const
{
    if (originalToNewCellMap_.found(originalCellId))
    {
        return originalToNewCellMap_[originalCellId];
    }
    
    return originalCellId;
}

void Foam::fvMeshTopoChangers::cellAddition::createOrdinaryBoundaryFaces(polyTopoChange& meshMod)
{
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    
    if (debugLevel_ >= 1)
    {
        Info<< "创建普通边界面: " << additionData_.ordinaryBoundaryFaces.size() << " 个面" << endl;
    }
    
    forAll(additionData_.ordinaryBoundaryFaces, i)
    {
        BoundaryFaceInfo& faceInfo = additionData_.ordinaryBoundaryFaces[i];
        
        const word& sourceBoundaryName = faceInfo.boundaryName;
        const word targetBoundaryName = getTargetBoundaryName(sourceBoundaryName);

        label patchID = -1;
        forAll(pbm, patchi)
        {
            if (pbm[patchi].name() == targetBoundaryName)
            {
                patchID = patchi;
                break;
            }
        }
        
        if (patchID == -1)
        {
            WarningInFunction
                << "找不到目标边界 '" << targetBoundaryName 
                << "' (源边界 '" << sourceBoundaryName << "')" << endl;
            continue;
        }
        
        // 构建面的点列表
        face f(faceInfo.pointIds.size());
        forAll(faceInfo.pointIds, j)
        {
            f[j] = getNewPointId(faceInfo.pointIds[j]);
        }
        
        label ownCellId = getNewCellId(faceInfo.ownCellId);
        
        const polyPatch& patch = pbm[patchID];
        label referenceFaceID = patch.start();

        label newFaceId = meshMod.addFace
        (
            f,
            ownCellId,
            -1,
            referenceFaceID,
            false,
            patchID
        );
        
        faceInfo.newId = newFaceId;
        
        if (debugLevel_ >= 2)
        {
            Info<< "  边界面 " << faceInfo.originalId << " -> " << newFaceId
                << " 边界: " << sourceBoundaryName << "->" << targetBoundaryName
                << " owner: " << faceInfo.ownCellId << "->" << ownCellId << endl;
        }
    }
}

void Foam::fvMeshTopoChangers::cellAddition::createInternalFaces(polyTopoChange& meshMod)
{
    if (debugLevel_ >= 1)
    {
        Info<< "创建内部面: " << additionData_.internalFaces.size() << " 个面" << endl;
    }
    
    forAll(additionData_.internalFaces, i)
    {
        InternalFaceInfo& faceInfo = additionData_.internalFaces[i];
        
        // 构建面的点列表
        face f(faceInfo.pointIds.size());
        forAll(faceInfo.pointIds, j)
        {
            f[j] = getNewPointId(faceInfo.pointIds[j]);
        }
        
        label ownCellId = getNewCellId(faceInfo.ownCellId);
        label neiCellId = getNewCellId(faceInfo.neiCellId);
        
        label newFaceId = meshMod.addFace
        (
            f,
            ownCellId,
            neiCellId,
            -1,
            false,
            -1
        );
        
        faceInfo.newId = newFaceId;
        
        if (debugLevel_ >= 2)
        {
            Info<< "  内部面 " << faceInfo.originalId << " -> " << newFaceId
                << " owner: " << faceInfo.ownCellId << "->" << ownCellId
                << " nei: " << faceInfo.neiCellId << "->" << neiCellId << endl;
        }
    }
}

void Foam::fvMeshTopoChangers::cellAddition::modifyExistingCouplingFaces(polyTopoChange& meshMod)
{
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const polyPatch& couplingPatch = pbm[couplingPatchID_];
    const labelList& faceOwner = mesh().faceOwner();
    const faceList& meshFaces = mesh().faces();
    
    if (debugLevel_ >= 1)
    {
        Info<< "修改现有耦合边界面为内部面..." << endl;
    }
    
    forAll(additionData_.couplingBoundaryInfo.faceIds, i)
    {
        const label patchFaceI = additionData_.couplingBoundaryInfo.faceIds[i];
        const bool ownerRemoved = additionData_.couplingBoundaryInfo.ownerRemoved[i];
        const label removedCellId = additionData_.couplingBoundaryInfo.removedCellIds[i];
        
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
            label newNeiCellId = getNewCellId(removedCellId);
            
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
                -1
            );
        }
    }
}

void Foam::fvMeshTopoChangers::cellAddition::createTemporaryCouplingFaces(polyTopoChange& meshMod)
{
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    
    if (debugLevel_ >= 1)
    {
        Info<< "创建临时耦合边界面: " << additionData_.newCouplingBoundaryFaces.size() << " 个面" << endl;
    }
    
    label temporaryPatchID = 0;
    const polyPatch& temporaryPatch = pbm[temporaryPatchID];
    label tempPatchCurrentSize = temporaryPatch.size();
    
    forAll(additionData_.newCouplingBoundaryFaces, i)
    {
        BoundaryFaceInfo& faceInfo = additionData_.newCouplingBoundaryFaces[i];
        
        face f(faceInfo.pointIds.size());
        forAll(faceInfo.pointIds, j)
        {
            f[j] = getNewPointId(faceInfo.pointIds[j]);
        }
        
        label ownCellId = getNewCellId(faceInfo.ownCellId);
        
        label newFaceId = meshMod.addFace
        (
            f,
            ownCellId,
            -1,
            -1,
            false,
            temporaryPatchID
        );
        
        faceInfo.newId = newFaceId;
        
        label expectedTempFaceId = tempPatchCurrentSize + i;
        temporaryBoundaryFaceIds_.append(expectedTempFaceId);
        
        if (debugLevel_ >= 2)
        {
            Info<< "  临时面 " << faceInfo.originalId << " -> " << newFaceId
                << " owner: " << faceInfo.ownCellId << "->" << ownCellId << endl;
        }
    }
}

void Foam::fvMeshTopoChangers::cellAddition::convertTemporaryToCouplingFaces(polyTopoChange& meshMod)
{
    if (debugLevel_ >= 1)
    {
        Info<< "将 " << temporaryBoundaryFaceIds_.size() << " 个临时面转换为耦合面" << endl;
    }
    
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const faceList& meshFaces = mesh().faces();
    const labelList& faceOwner = mesh().faceOwner();
    
    label temporaryPatchID = 0;
    const polyPatch& temporaryPatch = pbm[temporaryPatchID];
    
    forAll(temporaryBoundaryFaceIds_, i)
    {
        label faceId = temporaryBoundaryFaceIds_[temporaryBoundaryFaceIds_.size() - 1 - i];
        
        if (faceId >= 0 && faceId < temporaryPatch.size())
        {
            label globalFaceId = temporaryPatch.start() + faceId;
            const face& f = meshFaces[globalFaceId];
            label owner = faceOwner[globalFaceId];
            
            face flippedFace(f);
            flippedFace.flip();
            
            meshMod.modifyFace
            (
                flippedFace,
                globalFaceId,
                owner,
                -1,
                false,
                couplingPatchID_
            );
            
            if (debugLevel_ >= 2)
            {
                Info<< "  转换面 " << i << ": 全局ID=" << globalFaceId 
                    << " owner=" << owner << endl;
            }
        }
    }
}

void Foam::fvMeshTopoChangers::cellAddition::updateCellZones(const polyTopoChangeMap& map)
{
    if (debugLevel_ >= 1)
    {
        Info<< "更新cellZones..." << endl;
    }
    
    const labelList& cellMap = map.cellMap();
    
    DynamicList<label> addedCells;
    forAll(cellMap, cellI)
    {
        if (cellMap[cellI] == -1)
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
    
    cellZoneList& cellZones = const_cast<cellZoneList&>(mesh().cellZones());
    
    label targetCellZoneID = -1;
    forAll(cellZones, zoneI)
    {
        if (cellZones[zoneI].name() == "gas")
        {
            targetCellZoneID = zoneI;
            break;
        }
    }
    
    if (targetCellZoneID == -1)
    {
        WarningInFunction << "找不到目标cellZone: gas" << endl;
        return;
    }

    cellZone& targetCellZone = cellZones[targetCellZoneID];
    
    labelList oldCells = targetCellZone;
    labelList newCells(oldCells.size() + addedCells.size());
    
    forAll(oldCells, i)
    {
        newCells[i] = oldCells[i];
    }
    
    forAll(addedCells, i)
    {
        newCells[oldCells.size() + i] = addedCells[i];
    }
    
    Foam::stableSort(newCells);
    targetCellZone.transfer(newCells);
    
    if (debugLevel_ >= 1)
    {
        Info<< "已将 " << addedCells.size() << " 个新单元添加到目标cellZone" << endl;
        Info<< "cellZone更新后单元数: " << targetCellZone.size() << endl;
    }
}

Foam::autoPtr<Foam::polyTopoChangeMap> 
Foam::fvMeshTopoChangers::cellAddition::generateMesh()
{
    originalToNewPointMap_.clear();
    originalToNewCellMap_.clear();
    temporaryBoundaryFaceIds_.clear();
    
    mesh().preChange();
    
    polyTopoChange meshMod(mesh());
    
    // 第一步：创建点和单元
    createPoints(meshMod);
    createCells(meshMod);
    createTemporaryCouplingFaces(meshMod);
    modifyExistingCouplingFaces(meshMod);
    createOrdinaryBoundaryFaces(meshMod);
    createInternalFaces(meshMod);
    
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh(), false);
    meshChangeStep_ = 1;
    mesh().topoChange(map);
    
    // 第二步：转换临时面为耦合面
    mesh().preChange();
    polyTopoChange meshMod2(mesh());
    convertTemporaryToCouplingFaces(meshMod2);
    
    autoPtr<polyTopoChangeMap> map2 = meshMod2.changeMesh(mesh(), false);
    meshChangeStep_ = 2;
    mesh().topoChange(map2);
    
    return map2;
}

void Foam::fvMeshTopoChangers::cellAddition::debugOutput() const
{
    if (debugLevel_ < 1) return;
    
    Info<< nl << "=== 单元添加调试信息 ===" << nl;
    Info<< "时间: " << mesh().time().name() << nl;
    Info<< "源区域: " << sourceRegionName_ << nl;
    Info<< "耦合边界: " << couplingPatchName_ << nl;
    Info<< "添加点数: " << additionData_.points.size() << nl;
    Info<< "添加单元数: " << additionData_.cells.size() << nl;
    Info<< "普通边界面数: " << additionData_.ordinaryBoundaryFaces.size() << nl;
    Info<< "新耦合边界面数: " << additionData_.newCouplingBoundaryFaces.size() << nl;
    Info<< "内部面数: " << additionData_.internalFaces.size() << nl;
    Info<< "点映射数: " << originalToNewPointMap_.size() << nl;
    Info<< "单元映射数: " << originalToNewCellMap_.size() << nl;
    Info<< "=========================" << nl << endl;
}

// void Foam::fvMeshTopoChangers::cellAddition::debugOutput() const
// {
//     if (debugLevel_ < 1) return;
    
//     Info<< nl << "=== 单元添加调试信息 ===" << nl;
//     Info<< "时间: " << mesh().time().name() << nl;
//     Info<< "源区域: " << sourceRegionName_ << nl;
//     Info<< "耦合边界: " << couplingPatchName_ << nl;
//     Info<< "添加点数: " << additionData_.points.size() << nl;
//     Info<< "添加单元数: " << additionData_.cells.size() << nl;
//     Info<< "普通边界面数: " << additionData_.ordinaryBoundaryFaces.size() << nl;
//     Info<< "新耦合边界面数: " << additionData_.newCouplingBoundaryFaces.size() << nl;
//     Info<< "内部面数: " << additionData_.internalFaces.size() << nl;
//     Info<< "点映射数: " << originalToNewPointMap_.size() << nl;
//     Info<< "单元映射数: " << originalToNewCellMap_.size() << nl;
//     Info<< "=========================" << nl << endl;

//     if (debugLevel_ >= 4)
//     {
//         const pointField& meshPoints = mesh().points();
//         const faceList& meshFaces = mesh().faces();
//         const labelList& faceOwner = mesh().faceOwner();
//         const labelList& faceNeighbour = mesh().faceNeighbour();
//         const polyBoundaryMesh& pbm = mesh().boundaryMesh();
        
//         Info<< nl << "=== CELLADDITION 详细调试信息 (级别4) ===" << nl;
        
//         // 1. 网格统计信息
//         Info<< nl << "当前网格统计信息:" << nl;
//         Info<< "  网格点数: " << mesh().nPoints() << nl;
//         Info<< "  网格单元数: " << mesh().nCells() << nl;
//         Info<< "  网格面数: " << mesh().nFaces() << nl;
//         Info<< "  内部面数: " << mesh().nInternalFaces() << nl;
//         Info<< "  边界patch数: " << pbm.size() << nl;
        
//         forAll(pbm, patchi)
//         {
//             Info<< "    patch" << patchi << ": " << pbm[patchi].name() 
//                 << " (类型=" << pbm[patchi].type() << ", 面数=" << pbm[patchi].size() << ")" << nl;
//         }
        
//         // 2. 添加的点信息
//         Info<< nl << "=== 添加的点详情 ===" << nl;
//         Info<< "添加点总数: " << additionData_.points.size() << nl;
//         forAll(additionData_.points, i)
//         {
//             const PointInfo& pt = additionData_.points[i];
//             Info<< "  点" << i << ": 原ID=" << pt.originalId 
//                 << " -> 新ID=" << pt.newId
//                 << " 坐标=" << pt.coord << nl;
//         }
        
//         // 3. 添加的单元信息
//         Info<< nl << "=== 添加的单元详情 ===" << nl;
//         Info<< "添加单元总数: " << additionData_.cells.size() << nl;
//         forAll(additionData_.cells, i)
//         {
//             const CellInfo& cell = additionData_.cells[i];
//             Info<< "  单元" << i << ": 原ID=" << cell.originalId 
//                 << " -> 新ID=" << cell.newId
//                 << " 中心=" << cell.centroid << nl;
//         }
        
//         // 4. 普通边界面详情
//         Info<< nl << "=== 普通边界面详情 ===" << nl;
//         Info<< "普通边界面总数: " << additionData_.ordinaryBoundaryFaces.size() << nl;
//         forAll(additionData_.ordinaryBoundaryFaces, i)
//         {
//             const BoundaryFaceInfo& face = additionData_.ordinaryBoundaryFaces[i];
            
//             Info<< "  普通边界面" << i << ":" << nl;
//             Info<< "    原ID=" << face.originalId << " -> 新ID=" << face.newId << nl;
//             Info<< "    源边界=" << face.boundaryName;
            
//             // 查找目标边界
//             word targetBoundary = getTargetBoundaryName(face.boundaryName);
//             if (targetBoundary != face.boundaryName)
//             {
//                 Info<< " -> 目标边界=" << targetBoundary;
//             }
//             Info<< nl;
            
//             Info<< "    原owner=" << face.ownCellId << " -> 新owner=" << getNewCellId(face.ownCellId) << nl;
//             Info<< "    原点列表=" << face.pointIds << nl;
            
//             // 输出映射后的点列表和坐标
//             labelList newPointIds(face.pointIds.size());
//             pointField facePointCoords(face.pointIds.size());
//             forAll(face.pointIds, j)
//             {
//                 newPointIds[j] = getNewPointId(face.pointIds[j]);
                
//                 // 查找点坐标
//                 bool foundCoord = false;
//                 forAll(additionData_.points, ptI)
//                 {
//                     if (additionData_.points[ptI].originalId == face.pointIds[j])
//                     {
//                         facePointCoords[j] = additionData_.points[ptI].coord;
//                         foundCoord = true;
//                         break;
//                     }
//                 }
//                 if (!foundCoord)
//                 {
//                     facePointCoords[j] = point::zero;
//                 }
//             }
            
//             Info<< "    新点列表=" << newPointIds << nl;
//             Info<< "    点坐标:" << nl;
//             forAll(facePointCoords, j)
//             {
//                 Info<< "      点" << face.pointIds[j] << "(" << newPointIds[j] 
//                     << "): " << facePointCoords[j] << nl;
//             }
            
//             // 计算面中心和法向量
//             if (facePointCoords.size() >= 3)
//             {
//                 point faceCentre = average(facePointCoords);
//                 vector faceNormal = vector::zero;
//                 for (label j = 1; j < facePointCoords.size() - 1; j++)
//                 {
//                     faceNormal += (facePointCoords[j] - facePointCoords[0]) ^ 
//                                   (facePointCoords[j+1] - facePointCoords[0]);
//                 }
//                 faceNormal /= 2.0;
                
//                 Info<< "    面中心=" << faceCentre << nl;
//                 Info<< "    面法向量=" << faceNormal << nl;
//                 Info<< "    面积=" << mag(faceNormal) << nl;
//             }
//         }
        
//         // 5. 新耦合边界面详情
//         Info<< nl << "=== 新耦合边界面详情 ===" << nl;
//         Info<< "新耦合边界面总数: " << additionData_.newCouplingBoundaryFaces.size() << nl;
//         forAll(additionData_.newCouplingBoundaryFaces, i)
//         {
//             const BoundaryFaceInfo& face = additionData_.newCouplingBoundaryFaces[i];
            
//             Info<< "  新耦合面" << i << ":" << nl;
//             Info<< "    原ID=" << face.originalId << " -> 新ID=" << face.newId << nl;
//             Info<< "    源边界=" << face.boundaryName << " -> 目标耦合边界=" << couplingPatchName_ << nl;
//             Info<< "    原owner=" << face.ownCellId << " -> 新owner=" << getNewCellId(face.ownCellId) << nl;
//             Info<< "    原点列表=" << face.pointIds << nl;
            
//             // 输出映射后的点列表和坐标
//             labelList newPointIds(face.pointIds.size());
//             pointField facePointCoords(face.pointIds.size());
//             forAll(face.pointIds, j)
//             {
//                 newPointIds[j] = getNewPointId(face.pointIds[j]);
                
//                 bool foundCoord = false;
//                 forAll(additionData_.points, ptI)
//                 {
//                     if (additionData_.points[ptI].originalId == face.pointIds[j])
//                     {
//                         facePointCoords[j] = additionData_.points[ptI].coord;
//                         foundCoord = true;
//                         break;
//                     }
//                 }
//                 if (!foundCoord)
//                 {
//                     facePointCoords[j] = point::zero;
//                 }
//             }
            
//             Info<< "    新点列表=" << newPointIds << nl;
//             Info<< "    点坐标:" << nl;
//             forAll(facePointCoords, j)
//             {
//                 Info<< "      点" << face.pointIds[j] << "(" << newPointIds[j] 
//                     << "): " << facePointCoords[j] << nl;
//             }
            
//             // 计算面中心和法向量
//             if (facePointCoords.size() >= 3)
//             {
//                 point faceCentre = average(facePointCoords);
//                 vector faceNormal = vector::zero;
//                 for (label j = 1; j < facePointCoords.size() - 1; j++)
//                 {
//                     faceNormal += (facePointCoords[j] - facePointCoords[0]) ^ 
//                                   (facePointCoords[j+1] - facePointCoords[0]);
//                 }
//                 faceNormal /= 2.0;
                
//                 Info<< "    面中心=" << faceCentre << nl;
//                 Info<< "    面法向量=" << faceNormal << nl;
//                 Info<< "    面积=" << mag(faceNormal) << nl;
//             }
            
//             // 输出临时patch信息
//             if (temporaryBoundaryFaceIds_.size() > i)
//             {
//                 Info<< "    临时patch面ID=" << temporaryBoundaryFaceIds_[i] << nl;
//             }
//         }
        
//         // 6. 内部面详情
//         Info<< nl << "=== 内部面详情 ===" << nl;
//         Info<< "内部面总数: " << additionData_.internalFaces.size() << nl;
//         forAll(additionData_.internalFaces, i)
//         {
//             const InternalFaceInfo& face = additionData_.internalFaces[i];
            
//             Info<< "  内部面" << i << ":" << nl;
//             Info<< "    原ID=" << face.originalId << " -> 新ID=" << face.newId << nl;
//             Info<< "    原owner=" << face.ownCellId << " -> 新owner=" << getNewCellId(face.ownCellId) << nl;
//             Info<< "    原neighbour=" << face.neiCellId << " -> 新neighbour=" << getNewCellId(face.neiCellId) << nl;
//             Info<< "    原点列表=" << face.pointIds << nl;
            
//             // 输出映射后的点列表和坐标
//             labelList newPointIds(face.pointIds.size());
//             pointField facePointCoords(face.pointIds.size());
//             forAll(face.pointIds, j)
//             {
//                 newPointIds[j] = getNewPointId(face.pointIds[j]);
                
//                 bool foundCoord = false;
//                 forAll(additionData_.points, ptI)
//                 {
//                     if (additionData_.points[ptI].originalId == face.pointIds[j])
//                     {
//                         facePointCoords[j] = additionData_.points[ptI].coord;
//                         foundCoord = true;
//                         break;
//                     }
//                 }
//                 if (!foundCoord)
//                 {
//                     facePointCoords[j] = point::zero;
//                 }
//             }
            
//             Info<< "    新点列表=" << newPointIds << nl;
//             Info<< "    点坐标:" << nl;
//             forAll(facePointCoords, j)
//             {
//                 Info<< "      点" << face.pointIds[j] << "(" << newPointIds[j] 
//                     << "): " << facePointCoords[j] << nl;
//             }
            
//             // 计算面中心和法向量
//             if (facePointCoords.size() >= 3)
//             {
//                 point faceCentre = average(facePointCoords);
//                 vector faceNormal = vector::zero;
//                 for (label j = 1; j < facePointCoords.size() - 1; j++)
//                 {
//                     faceNormal += (facePointCoords[j] - facePointCoords[0]) ^ 
//                                   (facePointCoords[j+1] - facePointCoords[0]);
//                 }
//                 faceNormal /= 2.0;
                
//                 Info<< "    面中心=" << faceCentre << nl;
//                 Info<< "    面法向量=" << faceNormal << nl;
//                 Info<< "    面积=" << mag(faceNormal) << nl;
//             }
//         }
        
//         // 7. 现有耦合边界修改信息
//         Info<< nl << "=== 现有耦合边界修改详情 ===" << nl;
//         Info<< "待修改的耦合边界面数: " << additionData_.couplingBoundaryInfo.faceIds.size() << nl;
        
//         forAll(additionData_.couplingBoundaryInfo.faceIds, i)
//         {
//             const label patchFaceI = additionData_.couplingBoundaryInfo.faceIds[i];
//             const bool ownerRemoved = additionData_.couplingBoundaryInfo.ownerRemoved[i];
//             const label removedCellId = additionData_.couplingBoundaryInfo.removedCellIds[i];
            
//             Info<< "  耦合面修改" << i << ":" << nl;
//             Info<< "    patch内面ID=" << patchFaceI << nl;
//             Info<< "    owner被移除=" << (ownerRemoved ? "是" : "否") << nl;
//             Info<< "    被移除的单元ID=" << removedCellId << nl;
            
//             if (ownerRemoved)
//             {
//                 label newNeiCellId = getNewCellId(removedCellId);
//                 Info<< "    新neighbour单元=" << newNeiCellId << nl;
//                 Info<< "    变化: 边界面 -> 内部面" << nl;
//             }
            
//             // 如果耦合patch存在，显示面的详细信息
//             if (couplingPatchID_ >= 0 && couplingPatchID_ < pbm.size())
//             {
//                 const polyPatch& couplingPatch = pbm[couplingPatchID_];
                
//                 if (patchFaceI >= 0 && patchFaceI < couplingPatch.size())
//                 {
//                     const label globalFaceId = couplingPatch.start() + patchFaceI;
//                     const face& f = meshFaces[globalFaceId];
//                     const label currentOwner = faceOwner[globalFaceId];
                    
//                     Info<< "    全局面ID=" << globalFaceId << nl;
//                     Info<< "    当前owner=" << currentOwner << nl;
//                     Info<< "    面点列表=" << f << nl;
                    
//                     // 显示面的几何信息
//                     point faceCentre = f.centre(meshPoints);
//                     vector faceArea = f.area(meshPoints);
                    
//                     Info<< "    面中心=" << faceCentre << nl;
//                     Info<< "    面法向量=" << faceArea << nl;
//                     Info<< "    面积=" << mag(faceArea) << nl;
                    
//                     // 显示点坐标
//                     Info<< "    面点坐标:" << nl;
//                     forAll(f, j)
//                     {
//                         Info<< "      点" << f[j] << ": " << meshPoints[f[j]] << nl;
//                     }
//                 }
//             }
//         }
        
//         // 8. 当前耦合patch详细信息
//         if (couplingPatchID_ >= 0 && couplingPatchID_ < pbm.size())
//         {
//             const polyPatch& couplingPatch = pbm[couplingPatchID_];
//             Info<< nl << "=== 当前耦合patch详细信息 ===" << nl;
//             Info<< "耦合patch名称: " << couplingPatch.name() << nl;
//             Info<< "耦合patch类型: " << couplingPatch.type() << nl;
//             Info<< "耦合patch面数: " << couplingPatch.size() << nl;
//             Info<< "耦合patch起始面ID: " << couplingPatch.start() << nl;
            
//             // 显示前几个和后几个耦合面的详细信息
//             const label nFacesToShow = min(5, couplingPatch.size());
//             Info<< "前" << nFacesToShow << "个耦合面详情:" << nl;
//             for (label i = 0; i < nFacesToShow; i++)
//             {
//                 const label globalFaceI = couplingPatch.start() + i;
//                 const face& f = meshFaces[globalFaceI];
//                 const label owner = faceOwner[globalFaceI];
                
//                 Info<< "  耦合面" << i << " (全局ID=" << globalFaceI << "):" << nl;
//                 Info<< "    owner单元=" << owner << nl;
//                 Info<< "    点列表=" << f << nl;
                
//                 point faceCentre = f.centre(meshPoints);
//                 vector faceArea = f.area(meshPoints);
                
//                 Info<< "    面中心=" << faceCentre << nl;
//                 Info<< "    面法向量=" << faceArea << nl;
//                 Info<< "    面积=" << mag(faceArea) << nl;
//             }
            
//             if (couplingPatch.size() > 10)
//             {
//                 Info<< "... 中间" << (couplingPatch.size() - 10) << "个面未显示 ..." << nl;
                
//                 const label startIdx = max(nFacesToShow, couplingPatch.size() - 5);
//                 Info<< "后5个耦合面详情:" << nl;
//                 for (label i = startIdx; i < couplingPatch.size(); i++)
//                 {
//                     const label globalFaceI = couplingPatch.start() + i;
//                     const face& f = meshFaces[globalFaceI];
//                     const label owner = faceOwner[globalFaceI];
                    
//                     Info<< "  耦合面" << i << " (全局ID=" << globalFaceI << "):" << nl;
//                     Info<< "    owner单元=" << owner << nl;
//                     Info<< "    点列表=" << f << nl;
                    
//                     point faceCentre = f.centre(meshPoints);
//                     vector faceArea = f.area(meshPoints);
                    
//                     Info<< "    面中心=" << faceCentre << nl;
//                     Info<< "    面法向量=" << faceArea << nl;
//                     Info<< "    面积=" << mag(faceArea) << nl;
//                 }
//             }
//         }
        
//         // 9. 点和单元映射表信息
//         Info<< nl << "=== 映射表详情 ===" << nl;
//         Info<< "点映射表大小: " << originalToNewPointMap_.size() << nl;
//         Info<< "单元映射表大小: " << originalToNewCellMap_.size() << nl;
        
//         if (originalToNewPointMap_.size() > 0)
//         {
//             Info<< "点映射详情 (前10个):" << nl;
//             label count = 0;
//             forAllConstIter(Map<label>, originalToNewPointMap_, iter)
//             {
//                 if (count >= 10) break;
//                 Info<< "  原点" << iter.key() << " -> 新点" << iter() << nl;
//                 count++;
//             }
//         }
        
//         if (originalToNewCellMap_.size() > 0)
//         {
//             Info<< "单元映射详情 (前10个):" << nl;
//             label count = 0;
//             forAllConstIter(Map<label>, originalToNewCellMap_, iter)
//             {
//                 if (count >= 10) break;
//                 Info<< "  原单元" << iter.key() << " -> 新单元" << iter() << nl;
//                 count++;
//             }
//         }
        
//         // 10. 临时面ID信息
//         Info<< nl << "=== 临时面ID信息 ===" << nl;
//         Info<< "临时面ID列表大小: " << temporaryBoundaryFaceIds_.size() << nl;
//         if (temporaryBoundaryFaceIds_.size() > 0)
//         {
//             Info<< "临时面ID列表: " << temporaryBoundaryFaceIds_ << nl;
//         }
        
//         Info<< nl << "=== CELLADDITION 详细调试完成 ===" << nl << endl;
//     }
// }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::cellAddition::cellAddition(fvMesh& mesh, const dictionary& dict)
:
    fvMeshTopoChanger(mesh),
    dict_(dict),
    debugLevel_(dict_.lookupOrDefault<label>("debugLevel", 1)),
    couplingPatchID_(-1),
    couplingPatchName_("gas_to_solid"),
    sourceRegionName_("solid"),
    timeIndex_(-1),
    changedSinceWrite_(false),
    meshChangeStep_(0)
{
    readDict();

    if (debugLevel_ >= 1)
    {
        Info<< "创建单元添加网格拓扑变化器:" << nl
            << "  源区域: " << sourceRegionName_ << nl
            << "  耦合边界: " << couplingPatchName_ << nl
            << "  调试级别: " << debugLevel_ << nl
            << "  边界映射数量: " << boundaryNameMapping_.size() << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::cellAddition::~cellAddition()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::cellAddition::update()
{
    if (timeIndex_ == mesh().time().timeIndex())
    {
        return false;
    }
    
    timeIndex_ = mesh().time().timeIndex();
    
    if (readAdditionData() && additionData_.isValid)
    {
        if (debugLevel_ >= 1)
        {
            Info<< "检测到单元添加数据，开始生成新网格单元..." << endl;
        }
        
        autoPtr<polyTopoChangeMap> map = generateMesh();
        changedSinceWrite_ = true;
        return true;
    }
    
    return false;
}

void Foam::fvMeshTopoChangers::cellAddition::topoChange(const polyTopoChangeMap& map)
{
    if (debugLevel_ >= 1)
    {
        Info<< "cellAddition::topoChange : 处理拓扑变化 (步骤=" << meshChangeStep_ << ")" << endl;
    }
    
    switch (meshChangeStep_)
    {
        case 1:
        {
            updateCellZones(map);
            break;
        }
        
        case 2:
        {
            debugOutput();

            // 清理数据
            additionData_.clear();
            originalToNewPointMap_.clear();
            originalToNewCellMap_.clear();
            temporaryBoundaryFaceIds_.clear();
            meshChangeStep_ = 0;
            
            break;
        }
    }
}

void Foam::fvMeshTopoChangers::cellAddition::mapMesh(const polyMeshMap& map)
{
    NotImplemented;
}

void Foam::fvMeshTopoChangers::cellAddition::distribute(const polyDistributionMap& map)
{
    // 分布式映射处理
}

bool Foam::fvMeshTopoChangers::cellAddition::write(const bool write) const
{
    return changedSinceWrite_;
}

// ************************************************************************* //