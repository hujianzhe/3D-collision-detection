# 3D-collision-detection
支持碰撞体: 射线,平面,AABB立方体,球,胶囊,三角形.生成碰撞体的参数对象至少需要成员 {x, y, z}  
intersect运算返回true为相交状态.  
cast碰撞运算返回null为不碰撞,返回对象 {distance, hit_point, hit_normal} 为碰撞,distance为0时说明碰撞体处于相交状态.  
