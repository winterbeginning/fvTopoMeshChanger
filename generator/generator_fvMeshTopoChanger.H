#ifndef generator_fvMeshTopoChanger_H
#define generator_fvMeshTopoChanger_H

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

// 从固相传来的数据结构
struct ReceivedPointInfo
{
    label originalId;
    point coord;
    label newId;  // 在当前域的新ID
    
    ReceivedPointInfo() : originalId(-1), coord(point::zero), newId(-1) {}
    ReceivedPointInfo(label id, const point& c) : originalId(id), coord(c), newId(-1) {}
};

struct ReceivedCellInfo
{
    label originalId;
    point centroid;      // 新增：单元中心坐标
    label newId;         // 创建后的新ID
    
    ReceivedCellInfo()
    :
        originalId(-1),
        centroid(point::zero),
        newId(-1)
    {}
    
    ReceivedCellInfo(label origId, const point& centre)
    :
        originalId(origId),
        centroid(centre),
        newId(-1)
    {}
};

struct ReceivedFaceInfo
{
    label originalId;
    List<label> pointIds;
    label ownCellId;
    label neiCellId;  // 对于内部面有效
    word boundaryName;
    label faceType;   // 0=普通边界面, 1=新top边界面, 2=内部面
    label newId;      // 在当前域的新ID
    
    ReceivedFaceInfo() : originalId(-1), ownCellId(-1), neiCellId(-1), faceType(0), newId(-1) {}
};

struct ReceivedTopBoundaryInfo
{
    List<label> faceIds;           // patch内局部面ID
    List<bool> ownerRemoved;       // owner是否被移除
    List<label> removedCellIds;    // 对应的单元ID
};

struct ReceivedMeltingData
{
    List<ReceivedPointInfo> points;
    List<ReceivedCellInfo> cells;
    List<ReceivedFaceInfo> ordinaryBoundaryFaces;
    List<ReceivedFaceInfo> newTopBoundaryFaces;
    List<ReceivedFaceInfo> internalFaces;
    ReceivedTopBoundaryInfo topBoundaryInfo;
    
    // 新增：ID映射关系
    Map<label> originalToNewPointMap;  // 原始点ID -> 新创建的点ID
    Map<label> originalToNewCellMap;   // 原始单元ID -> 新创建的单元ID

    bool isValid;
    word sourceTime;
    
    ReceivedMeltingData() : isValid(false) {}
    
    void clear()
    {
        points.clear();
        cells.clear();
        ordinaryBoundaryFaces.clear();
        newTopBoundaryFaces.clear();
        internalFaces.clear();
        topBoundaryInfo.faceIds.clear();
        topBoundaryInfo.ownerRemoved.clear();
        topBoundaryInfo.removedCellIds.clear();
        originalToNewPointMap.clear();  // 新增
        originalToNewCellMap.clear();   // 新增        
        isValid = false;
        sourceTime = word::null;
    }
    
    bool readFromGlobalFile(const fileName& dataFile);
    void debugPrintFileContent(const fileName& dataFile);
};

/*---------------------------------------------------------------------------*\
                        Class generator Declaration
\*---------------------------------------------------------------------------*/

class generator
:
    public fvMeshTopoChanger
{
    // Private Data

        //- 配置字典
        dictionary dict_;
        
        label debugLevel_; 

            // ID映射表
        Map<label> originalToNewPointMap_;
        Map<label> originalToNewCellMap_;
        //- 固相区域名称
        word solidRegionName_;
        
        //- 耦合边界补丁ID
        label couplingPatchID_;
        
        //- 耦合边界名称
        word couplingPatchName_;

        //- 边界名称映射表：固相边界名 -> 气相边界名
        HashTable<word, word> boundaryNameMapping_; 
        
        //- 接收到的融化数据
        ReceivedMeltingData meltingData_;
        
        //- 自上次写入以来是否发生了拓扑变化
        bool changedSinceWrite_;
        
        //- 当前时间步索引
        label timeIndex_;

        //- 临时边界面ID列表（用于两步网格变化）
        DynamicList<label> temporaryBoundaryFaceIds_;

        label meshChangeStep_;

    // Private Member Functions

        //- 从字典中读取参数
        void readDict();

        //- 读取边界名称映射
        void readBoundaryMapping();
        
        //- 找到耦合边界ID
        void findCouplingPatch();

        //- 根据固相边界名获取对应的气相边界名
        word getGasBoundaryName(const word& solidBoundaryName) const;
        
        //- 尝试读取固相传来的融化数据
        bool readMeltingData();
        
        //- 根据接收到的数据生成新网格
        autoPtr<polyTopoChangeMap> generateMesh();
        
        //- 创建点并返回新的点ID
        void createPoints(polyTopoChange& meshMod);
        
        //- 创建单元并返回新的单元ID
        void createCells(polyTopoChange& meshMod);
        
        //- 创建普通边界面
        void createOrdinaryBoundaryFaces(polyTopoChange& meshMod);
        
        //- 创建新的top边界面（耦合面）
        //void createNewTopBoundaryFaces(polyTopoChange& meshMod);
        
        //- 创建内部面
        void createInternalFaces(polyTopoChange& meshMod);
        
        //- 修改现有耦合边界面为内部面
        void modifyExistingCouplingFaces(polyTopoChange& meshMod);

        label getNewPointId(label originalPointId) const;
        label getNewCellId(label originalCellId) const;

        void updateCellZones(const polyTopoChangeMap& map);


        //- 创建新的top边界面作为临时边界面
        void createNewTopBoundaryFacesAsTemporary(polyTopoChange& meshMod);
        
        //- 将临时边界面转换为耦合边界面
        void convertTemporaryFacesToCouplingFaces(polyTopoChange& meshMod);
        
        void customFieldMapping(const polyTopoChangeMap& map);
        
        //- 映射标量场
        void mapScalarFields(const labelList& newCells);
        
        //- 映射矢量场
        void mapVectorFields(const labelList& newCells);
        
        //- 从相邻单元插值标量
        scalar interpolateScalarFromNeighbors
        (
            const volScalarField& field,
            const label newCelli
        ) const;
        
        //- 从相邻单元插值矢量
        vector interpolateVectorFromNeighbors
        (
            const volVectorField& field,
            const label newCelli
        ) const;
        
        //- 获取场的默认值
        scalar getDefaultValue(const word& fieldName) const;
        
        //- 检查值是否合理
        bool isReasonableValue(const word& fieldName, const scalar value) const;

        //- 验证映射结果
        void validateMappingResults(const labelList& newCells);

            //- 基于单元中心坐标选择最佳母单元
        label selectBestMasterCell(const point& newCellCentre) const;

        // 调试输出函数
        void debugOutput(label level) const;
        
        //- 禁止默认拷贝构造
        generator(const generator&) = delete;
        
        //- 禁止赋值操作符
        void operator=(const generator&) = delete;

public:

    //- 运行时类型信息
    TypeName("generator");

    // Constructors

        //- 从网格和字典构造
        generator(fvMesh& mesh, const dictionary& dict);

    //- 析构函数
    virtual ~generator();

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