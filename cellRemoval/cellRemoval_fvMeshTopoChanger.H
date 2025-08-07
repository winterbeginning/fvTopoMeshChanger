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
    Foam::fvMeshTopoChangers::cellRemoval

Description
    基于温度阈值的单元移除网格拓扑变化器。
    移除温度超过指定阈值的单元，并在耦合界面处建立点映射关系。

SourceFiles
    cellRemoval_fvMeshTopoChanger.C

\*---------------------------------------------------------------------------*/

#ifndef cellRemoval_fvMeshTopoChanger_H
#define cellRemoval_fvMeshTopoChanger_H

#include "fvMeshTopoChanger.H"
#include "DynamicList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                        Class cellRemoval Declaration
\*---------------------------------------------------------------------------*/

class cellRemoval
:
    public fvMeshTopoChanger
{
public:

    // Public data structures

        //- 点信息结构
        struct PointInfo
        {
            label id;
            point coord;
            
            PointInfo() : id(-1), coord(point::zero) {}
            PointInfo(label pointId, const point& pt)
            : id(pointId), coord(pt) {}
        };

        //- 单元信息结构
        struct CellInfo
        {
            label cellId;
            point centroid;
            
            CellInfo() : cellId(-1), centroid(point::zero) {}
            CellInfo(label id, const point& center)
            : cellId(id), centroid(center) {}
        };

        //- 面信息基类
        struct FaceInfo
        {
            label faceId;
            labelList pointIds;
            label ownCellId;
            word faceType;
            
            FaceInfo() : faceId(-1), pointIds(), ownCellId(-1), faceType("face") {}
            FaceInfo(label id, const labelList& points, label owner, const word& type = "face")
            : faceId(id), pointIds(points), ownCellId(owner), faceType(type) {}
        };

        //- 边界面信息
        struct BoundaryFaceInfo : public FaceInfo
        {
            word boundaryName;
            
            BoundaryFaceInfo() : FaceInfo(), boundaryName("boundary") {}
            BoundaryFaceInfo(label id, const labelList& points, label owner, const word& boundary)
            : FaceInfo(id, points, owner, "boundary"), boundaryName(boundary) {}
        };

        //- 内部面信息
        struct InternalFaceInfo : public FaceInfo
        {
            label neiCellId;
            
            InternalFaceInfo() : FaceInfo(), neiCellId(-1) {}
            InternalFaceInfo(label id, const labelList& points, label owner, label neighbor)
            : FaceInfo(id, points, owner, "internal"), neiCellId(neighbor) {}
        };

        //- 数据输出结构
        struct RemovalDataOutput
        {
            // 调试级别
            label debugLevel_ = 1;
            
            // 移除的单元和点
            DynamicList<CellInfo> removedCells;
            DynamicList<PointInfo> removedCellPoints;
            
            // 耦合边界信息
            DynamicList<label> couplingBoundaryPoints;
            DynamicList<label> couplingBoundaryFaceIds;
            DynamicList<bool> couplingFaceOwnerRemoved;
            DynamicList<label> couplingFaceRemovedCellIds;
            
            // 面信息
            DynamicList<BoundaryFaceInfo> ordinaryBoundaryFaces;
            DynamicList<BoundaryFaceInfo> newCouplingBoundaryFaces;
            DynamicList<InternalFaceInfo> internalFaces;
            
            // 耦合点映射
            Map<label> couplingPointMapping;
            
            void clear();
            void applyCouplingPointMapping();
            void writeToFile(const polyMesh& mesh) const;
        };

private:

    // Private Data

        //- 字典
        const dictionary& dict_;

        //- 调试级别
        label debugLevel_;

        //- 融化温度阈值
        scalar temperatureThreshold_;

        //- 变更间隔
        label changeInterval_;

        //- 耦合边界ID
        label couplingPatchID_;

        //- 温度场名称
        word temperatureFieldName_;

        //- 是否已修改
        bool changedSinceWrite_;

        //- 时间索引
        label timeIndex_;

        //- 移除的单元列表
        labelList removedCells_;

    // Private Member Functions

        //- 读取字典
        void readDict();

        //- 查找耦合边界
        void findCouplingPatch();

        //- 选择要移除的单元
        labelList selectCellsToRemove() const;

        //- 建立耦合点映射
        void buildCouplingPointMapping(RemovalDataOutput& outputData) const;

        //- 执行网格变更
        autoPtr<polyTopoChangeMap> changeMesh(const labelList& cellsToRemove);

        //- 调试输出
        void debugOutput(const RemovalDataOutput& outputData) const;

public:

    //- Runtime type information
    TypeName("cellRemoval");

    // Constructors

        //- Construct from fvMesh and dictionary
        cellRemoval(fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        cellRemoval(const cellRemoval&) = delete;

    //- Destructor
    virtual ~cellRemoval();

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
        void operator=(const cellRemoval&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //