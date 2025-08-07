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
            Pout<< "读取边界名称映射:" << nl;
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
                    Pout<< "  " << sourceBoundaryName << " -> " << targetBoundaryName << nl;
                }
            }
        }
        
        if (debugLevel_ >= 1)
        {
            Pout<< "共读取 " << boundaryNameMapping_.size() << " 个边界映射" << endl;
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
    fileName globalDataDir;
    
    if (Pstream::parRun())
    {
        globalDataDir = runTime.rootPath()/runTime.globalCaseName()/"cellRemovalData"/("processor" + Foam::name(Pstream::myProcNo()));
    }
    else
    {
        globalDataDir = runTime.rootPath()/runTime.globalCaseName()/"cellRemovalData";
    }

    fileName dataFile = globalDataDir/("removalData_" + runTime.name() + ".dat");

    additionData_.clear();
    
    if (additionData_.readFromFile(dataFile))
    {
        additionData_.sourceTime = runTime.name();
        
        if (debugLevel_ >= 1)
        {
            Pout<< "成功读取单元添加数据: " << dataFile << nl
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
        Pout<< "创建点: " << additionData_.points.size() << " 个新点" << endl;
    }
    
    forAll(additionData_.points, i)
    {
        PointInfo& ptInfo = additionData_.points[i];
        
        label newPointId = meshMod.addPoint(ptInfo.coord, -1, true);
        
        originalToNewPointMap_.insert(ptInfo.originalId, newPointId);
        ptInfo.newId = newPointId;
        
        if (debugLevel_ >= 3)
        {
            Pout<< "  点 " << ptInfo.originalId << " -> " << newPointId 
                << " 坐标: " << ptInfo.coord << endl;
        }
    }
    
    if (debugLevel_ >= 1)
    {
        Pout<< "点创建完成，建立了 " << originalToNewPointMap_.size() << " 个点映射" << endl;
    }
}

void Foam::fvMeshTopoChangers::cellAddition::createCells(polyTopoChange& meshMod)
{
    if (debugLevel_ >= 1)
    {
        Pout<< "创建单元: " << additionData_.cells.size() << " 个新单元" << endl;
    }
    
    forAll(additionData_.cells, i)
    {
        CellInfo& cellInfo = additionData_.cells[i];

        label newCellId = meshMod.addCell(-1);
        
        originalToNewCellMap_.insert(cellInfo.originalId, newCellId);
        cellInfo.newId = newCellId;
        
        if (debugLevel_ >= 3)
        {
            Pout<< "  单元 " << cellInfo.originalId << " -> " << newCellId << endl;
        }
    }
    
    if (debugLevel_ >= 1)
    {
        Pout<< "单元创建完成，建立了 " << originalToNewCellMap_.size() << " 个单元映射" << endl;
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
        Pout<< "创建普通边界面: " << additionData_.ordinaryBoundaryFaces.size() << " 个面" << endl;
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
            Pout<< "  边界面 " << faceInfo.originalId << " -> " << newFaceId
                << " 边界: " << sourceBoundaryName << "->" << targetBoundaryName
                << " owner: " << faceInfo.ownCellId << "->" << ownCellId << endl;
        }
    }
}

void Foam::fvMeshTopoChangers::cellAddition::createInternalFaces(polyTopoChange& meshMod)
{
    if (debugLevel_ >= 1)
    {
        Pout<< "创建内部面: " << additionData_.internalFaces.size() << " 个面" << endl;
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
            Pout<< "  内部面 " << faceInfo.originalId << " -> " << newFaceId
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
        Pout<< "修改现有耦合边界面为内部面..." << endl;
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
                Pout<< "  修改耦合面 " << patchFaceI << " (全局ID: " << globalFaceId << ")"
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
        Pout<< "创建临时耦合边界面: " << additionData_.newCouplingBoundaryFaces.size() << " 个面" << endl;
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
            Pout<< "  临时面 " << faceInfo.originalId << " -> " << newFaceId
                << " owner: " << faceInfo.ownCellId << "->" << ownCellId << endl;
        }
    }
}

void Foam::fvMeshTopoChangers::cellAddition::convertTemporaryToCouplingFaces(polyTopoChange& meshMod)
{
    if (debugLevel_ >= 1)
    {
        Pout<< "将 " << temporaryBoundaryFaceIds_.size() << " 个临时面转换为耦合面" << endl;
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
                Pout<< "  转换面 " << i << ": 全局ID=" << globalFaceId 
                    << " owner=" << owner << endl;
            }
        }
    }
}

void Foam::fvMeshTopoChangers::cellAddition::updateCellZones(const polyTopoChangeMap& map)
{
    if (debugLevel_ >= 1)
    {
        Pout<< "更新cellZones..." << endl;
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
            Pout<< "没有新创建的单元，跳过cellZone更新" << endl;
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
        Pout<< "已将 " << addedCells.size() << " 个新单元添加到目标cellZone" << endl;
        Pout<< "cellZone更新后单元数: " << targetCellZone.size() << endl;
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
    
    Pout<< nl << "=== 单元添加调试信息 ===" << nl;
    Pout<< "时间: " << mesh().time().name() << nl;
    Pout<< "源区域: " << sourceRegionName_ << nl;
    Pout<< "耦合边界: " << couplingPatchName_ << nl;
    Pout<< "添加点数: " << additionData_.points.size() << nl;
    Pout<< "添加单元数: " << additionData_.cells.size() << nl;
    Pout<< "普通边界面数: " << additionData_.ordinaryBoundaryFaces.size() << nl;
    Pout<< "新耦合边界面数: " << additionData_.newCouplingBoundaryFaces.size() << nl;
    Pout<< "内部面数: " << additionData_.internalFaces.size() << nl;
    Pout<< "点映射数: " << originalToNewPointMap_.size() << nl;
    Pout<< "单元映射数: " << originalToNewCellMap_.size() << nl;
    Pout<< "=========================" << nl << endl;
}

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
        Pout<< "创建单元添加网格拓扑变化器:" << nl
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
            Pout<< "检测到单元添加数据，开始生成新网格单元..." << endl;
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
        Pout<< "cellAddition::topoChange : 处理拓扑变化 (步骤=" << meshChangeStep_ << ")" << endl;
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