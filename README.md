# 3D-collision-detection
支持碰撞体:射线,平面,AABB立方体,球,胶囊,三角形(部分支持).
碰撞运算返回对象 {distance, hit_point, hit_normal}, distance为0时说明碰撞体处于相交状态.
内部向量运算可自定义,但推荐使用three.js的Vector3类,自定义的需要实现基本向量操作来支持内部运算,接口形式要与three.js的一致.(待更新)
