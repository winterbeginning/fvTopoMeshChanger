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

Class
    Foam::fvMeshTopoChangers::cellAddition

Description
    基于外部数据的单元添加网格拓扑变化器。
    根据从其他区域传来的数据添加新的单元、点和面。

SourceFiles
    cellAddition_fvMeshTopoChanger.C

\*---------------------------------------------------------------------------*/

#ifndef cellAddition_fvMeshTopoChanger_H
#define cellAddition_fvMeshTopoChanger_H

#include "fvMeshTopoChanger.H"
#include "PackedBoolList.H"
#include "polyTopoChange.H"
#include "volFields.H"
#include "Time.H"
#include "IFstream.H"
#include "fileName.H"
#include "HashTable.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                        Class cellAddition Declaration
\*---------------------------------------------------------------------------*/

class cellAddition
:
    public fvMeshTopoChanger
{
public:

    // Public data structures

        //- 接收的点信息
        struct PointInfo
        {
            label originalId;
            point coord;
            label newId;
            
            PointInfo() : originalId(-1), coord(point::zero), newId(-1) {}
            PointInfo(label id, const point& c) : originalId(id), coord(c), newId(-1) {}
        };

        //- 接收的单元信息
        struct CellInfo
        {
            label originalId;
            point centroid;
            label newId;
            
            CellInfo() : originalId(-1), centroid(point::zero), newId(-1) {}
            CellInfo(label origId, const point& centre) : originalId(origId), centroid(centre), newId(-1) {}
        };

        //- 面信息基类
        struct FaceInfo
        {
            label originalId;
            labelList pointIds;
            label ownCellId;
            label neiCellId;
            word boundaryName;
            label faceType;
            label newId;
            
            FaceInfo() : originalId(-1), ownCellId(-1), neiCellId(-1), faceType(0), newId(-1) {}
        };

        //- 边界面信息
        struct BoundaryFaceInfo : public FaceInfo
        {
            BoundaryFaceInfo() : FaceInfo() {}
            BoundaryFaceInfo(label id, const labelList& points, label owner, const word& boundary)
            {
                originalId = id;
                pointIds = points;
                ownCellId = owner;
                boundaryName = boundary;
                faceType = 0;
                newId = -1;
            }
        };

        //- 内部面信息
        struct InternalFaceInfo : public FaceInfo
        {
            InternalFaceInfo() : FaceInfo() { faceType = 2; }
            InternalFaceInfo(label id, const labelList& points, label owner, label neighbor)
            {
                originalId = id;
                pointIds = points;
                ownCellId = owner;
                neiCellId = neighbor;
                faceType = 2;
                newId = -1;
            }
        };

        //- 耦合边界信息
        struct CouplingBoundaryInfo
        {
            DynamicList<label> faceIds;
            DynamicList<bool> ownerRemoved;
            DynamicList<label> removedCellIds;
            
            void clear()
            {
                faceIds.clear();
                ownerRemoved.clear();
                removedCellIds.clear();
            }
        };

        //- 数据接收结构
        struct AdditionDataInput
        {
            DynamicList<PointInfo> points;
            DynamicList<CellInfo> cells;
            DynamicList<BoundaryFaceInfo> ordinaryBoundaryFaces;
            DynamicList<BoundaryFaceInfo> newCouplingBoundaryFaces;
            DynamicList<InternalFaceInfo> internalFaces;
            CouplingBoundaryInfo couplingBoundaryInfo;
            
            bool isValid;
            word sourceTime;
            
            AdditionDataInput() : isValid(false) {}
            
            void clear();
            bool readFromFile(const fileName& dataFile);
        };

private:

    // Private Data

        //- 字典
        const dictionary& dict_;

        //- 调试级别
        label debugLevel_;

        //- 耦合边界ID
        label couplingPatchID_;

        //- 耦合边界名称
        word couplingPatchName_;

        //- 源区域名称
        word sourceRegionName_;

        //- 边界名称映射表
        HashTable<word, word> boundaryNameMapping_;

        //- 接收的数据
        AdditionDataInput additionData_;

        //- ID映射表
        Map<label> originalToNewPointMap_;
        Map<label> originalToNewCellMap_;

        //- 时间索引
        label timeIndex_;

        //- 是否已修改
        bool changedSinceWrite_;

        //- 网格变化步骤
        label meshChangeStep_;

        //- 临时边界面ID列表
        DynamicList<label> temporaryBoundaryFaceIds_;

    // Private Member Functions

        //- 读取字典
        void readDict();

        //- 查找耦合边界
        void findCouplingPatch();

        //- 读取边界名称映射
        void readBoundaryMapping();

        //- 获取目标边界名称
        word getTargetBoundaryName(const word& sourceBoundaryName) const;

        //- 读取数据文件
        bool readAdditionData();

        //- 生成网格
        autoPtr<polyTopoChangeMap> generateMesh();

        //- 创建点
        void createPoints(polyTopoChange& meshMod);

        //- 创建单元
        void createCells(polyTopoChange& meshMod);

        //- 创建普通边界面
        void createOrdinaryBoundaryFaces(polyTopoChange& meshMod);

        //- 创建内部面
        void createInternalFaces(polyTopoChange& meshMod);

        //- 修改现有耦合面
        void modifyExistingCouplingFaces(polyTopoChange& meshMod);

        //- 创建临时边界面
        void createTemporaryCouplingFaces(polyTopoChange& meshMod);

        //- 转换临时面为耦合面
        void convertTemporaryToCouplingFaces(polyTopoChange& meshMod);

        //- 获取新点ID
        label getNewPointId(label originalPointId) const;

        //- 获取新单元ID
        label getNewCellId(label originalCellId) const;

        //- 更新单元区域
        void updateCellZones(const polyTopoChangeMap& map);

        //- 调试输出
        void debugOutput() const;

public:

    //- Runtime type information
    TypeName("cellAddition");

    // Constructors

        //- Construct from fvMesh and dictionary
        cellAddition(fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        cellAddition(const cellAddition&) = delete;

    //- Destructor
    virtual ~cellAddition();

    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap& map);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap& map);

        //- Write using given format, version and compression
        virtual bool write(const bool write) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellAddition&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //