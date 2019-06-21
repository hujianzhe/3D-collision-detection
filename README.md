# 3D-collision-detection
支持碰撞体: 射线,平面,AABB立方体,球,胶囊,三角形(部分支持).生成碰撞体的参数对象至少需要成员 {x, y, z}  
碰撞运算返回对象 {distance, hit_point, hit_normal}, distance为0时说明碰撞体处于相交状态.   
