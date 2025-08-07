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
    Foam::fvMeshTopoChangers::melter

Description
    温度驱动的固相融化网格拓扑变化器
    
    移除温度超过融化温度的固相单元，并将相应表面转换为热通量边界

SourceFiles
    melter_fvMeshTopoChanger.C

\*---------------------------------------------------------------------------*/

#ifndef melter_fvMeshTopoChanger_H
#define melter_fvMeshTopoChanger_H

#include "fvMeshTopoChanger.H"
#include "PackedBoolList.H"
#include "polyTopoChange.H"
#include "volFields.H"
#include "Time.H"          // 新增
#include "OFstream.H"      // 新增
#include "fileName.H"      // 新增

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

// 完善数据结构
struct PointInfo
{
    label id;
    point coord;
    
    // 构造函数
    PointInfo() : id(-1), coord(point::zero) {}
    PointInfo(label i, const point& c) : id(i), coord(c) {}
};

// 在melter_fvMeshTopoChanger.H中修改CellInfo结构
struct CellInfo
{
    label cellId;
    point centroid;  // 新增：单元中心坐标
    
    CellInfo() : cellId(-1), centroid(point::zero) {}
    CellInfo(label id, const point& centre) :cellId(id), centroid(centre){}
};

// 索引映射结构
struct IndexMapping
{
    Map<label> oldToNewPointMap;     // 旧点ID -> 新点ID映射
    Map<label> oldToNewCellMap;      // 旧单元ID -> 新单元ID映射  
    Map<label> oldToNewFaceMap;      // 旧面ID -> 新面ID映射
    
    labelList newPointList;          // 新点ID列表
    labelList newCellList;           // 新单元ID列表
    labelList newFaceList;           // 新面ID列表
    
    void clear()
    {
        oldToNewPointMap.clear();
        oldToNewCellMap.clear();
        oldToNewFaceMap.clear();
        newPointList.clear();
        newCellList.clear();
        newFaceList.clear();
    }
};

// 修改面信息结构，添加类型标识
struct BoundaryFaceInfo
{
    label faceId;
    List<label> pointIds;
    label ownCellId;
    word boundaryName;
    label faceType;  // 新增：面类型标识 (0=普通边界面)
    
    // 构造函数
    BoundaryFaceInfo() : faceId(-1), ownCellId(-1), faceType(0) {}
    BoundaryFaceInfo(label id, const List<label>& pts, label own, const word& name, label type = 0) 
        : faceId(id), pointIds(pts), ownCellId(own), boundaryName(name), faceType(type) {}
};

struct NewSolidToGasBoundaryFaceInfo
{
    label faceId;
    List<label> pointIds;
    label ownCellId;
    word boundaryName;
    label faceType;  // 面类型标识 (1=新solid_to_gas边界面)
    
    // 构造函数
    NewSolidToGasBoundaryFaceInfo() : faceId(-1), ownCellId(-1), boundaryName("solid_to_gas"), faceType(1) {}
    NewSolidToGasBoundaryFaceInfo(label id, const List<label>& pts, label own, label type = 1) 
        : faceId(id), pointIds(pts), ownCellId(own), boundaryName("solid_to_gas"), faceType(type) {}
};

struct InternalFaceInfo
{
    label faceId;
    List<label> pointIds;
    label ownCellId;
    label neiCellId;
    label faceType;  // 新增：面类型标识 (2=内部面)
    
    // 构造函数
    InternalFaceInfo() : faceId(-1), ownCellId(-1), neiCellId(-1), faceType(2) {}
    InternalFaceInfo(label id, const List<label>& pts, label own, label nei, label type = 2) 
        : faceId(id), pointIds(pts), ownCellId(own), neiCellId(nei), faceType(type) {}
};

// 修改输出数据结构，添加重排功能
struct MeltingDataOutput
{
    label   debugLevel_;  // 调试级别
    // 原有数据
    List<PointInfo> removedCellPoints;
    List<CellInfo> removedCells;
    List<label> solidToGasBoundaryPoints;
    
    List<bool> solidToGasFaceOwnerRemoved;  // 修改：从topFaceOwnerRemoved
    List<label> solidToGasBoundaryFaceIds;  // 修改：从topBoundaryFaceIds
    List<label> solidToGasFaceRemovedCellIds;  // 修改：从topFaceRemovedCellIds

    List<BoundaryFaceInfo> ordinaryBoundaryFaces;
    List<NewSolidToGasBoundaryFaceInfo> newSolidToGasBoundaryFaces;
    List<InternalFaceInfo> internalFaces;
    
    // 新增：点索引映射（用于耦合界面）
    Map<label> solidToGasPointMapping;  // 固相界面点 -> 气相界面对应点

    // 新增：重排相关数据
    IndexMapping indexMapping;
    bool isReordered;  // 标记是否已重排
    
    void clear()
    {
         // 清空原有数据
        removedCellPoints.clear();
        removedCells.clear();
        solidToGasBoundaryPoints.clear();  // 修改
        solidToGasFaceOwnerRemoved.clear();  // 修改
        solidToGasBoundaryFaceIds.clear();  // 修改
        solidToGasFaceRemovedCellIds.clear();  // 修改
        ordinaryBoundaryFaces.clear();
        newSolidToGasBoundaryFaces.clear();  // 修改
        internalFaces.clear();
        
        // 清空映射数据
        solidToGasPointMapping.clear();  // 新增
        
        // 清空重排数据
        indexMapping.clear();
        isReordered = false;
    }
    
    void applyCouplingPointMapping();  // 新增：重排索引函数
    void writeToGlobalFile(const polyMesh& mesh) const;  // 修改：添加mesh参数
};


/*---------------------------------------------------------------------------*\
                        Class melter Declaration
\*---------------------------------------------------------------------------*/

class melter
:
    public fvMeshTopoChanger
{
    // Private Data

        //- 配置字典
        dictionary dict_;

        //- 调试级别
        label debugLevel_;
        
        //- 融化温度阈值
        scalar meltingTemperature_;
        
        //- 每隔多少时间步执行一次
        label changeInterval_;
        
        //- 热通量边界补丁ID
        label couplingPatchID_; 
        
        //- 温度场名称
        word temperatureFieldName_;
        
        //- 自上次写入以来是否发生了拓扑变化
        bool changedSinceWrite_;
        
        //- 当前时间步索引
        label timeIndex_;
        
        //- 上一次移除的单元
        labelList removedCells_;

    // Private Member Functions

        void collectAndOutputMeltingData(const labelList& cellsToRemove);
    
        //- 从字典中读取参数
        void readDict();
        
        //- 找到热通量边界ID
        void findSolidToGasPatch(); 
        
        //- 根据温度场选择要移除的单元
        labelList selectCellsToRemove() const;
        
        //- 执行网格变更
        autoPtr<polyTopoChangeMap> changeMesh(const labelList& cellsToRemove);
        
        // 新增：获取耦合界面点映射
        void buildCouplingPointMapping(MeltingDataOutput& outputData) const;

        void debugOutput(const MeltingDataOutput& outputData, label level) const;

        //- 禁止默认拷贝构造
        melter(const melter&) = delete;
        
        //- 禁止赋值操作符
        void operator=(const melter&) = delete;


public:

    //- 运行时类型信息
    TypeName("melter");


    // Constructors

        //- 从网格和字典构造
        melter(fvMesh& mesh, const dictionary& dict);


    //- 析构函数
    virtual ~melter();


    // Member Functions

        //- 更新网格拓扑（实现父类纯虚函数）
        virtual bool update();

        //- 响应拓扑变化（实现父类纯虚函数）
        virtual void topoChange(const polyTopoChangeMap&);

        //- 从另一个网格映射（实现父类纯虚函数）
        virtual void mapMesh(const polyMeshMap&);

        //- 响应分布式映射（实现父类纯虚函数）
        virtual void distribute(const polyDistributionMap&);

        //- 写出数据（实现父类纯虚函数）
        virtual bool write(const bool write = true) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //