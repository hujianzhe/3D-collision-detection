var CCT = {
	COLLISION_BODY_RAY : 0,
	COLLISION_BODY_AABB : 1,
	COLLISION_BODY_SPHERE : 2,
	COLLISION_BODY_CAPSULE : 3,
	COLLISION_BODY_PLANE : 4,
	COLLISION_BODY_TRIANGLES_PLANE : 5,
	
	EPSILON : 1E-5,
	Vector3 : function(x, y, z) {
		this.x = x || 0;
		this.y = y || 0;
		this.z = z || 0;
	},
	
	Ray : function (pos) {
		var o = {};
		o.type = CCT.COLLISION_BODY_RAY;
		if (pos)
			o.pos = new CCT.Vector3(pos.x, pos.y, pos.z);
		return o;
	},
	
	AABB : function (pos, half) {
		var o = {};
		o.type = CCT.COLLISION_BODY_AABB;
		if (pos)
			o.pos = new CCT.Vector3(pos.x, pos.y, pos.z);
		if (half)
			o.half = new CCT.Vector3(half.x, half.y, half.z);
		return o;
	},
	
	Sphere : function (pos, radius) {
		var o = {};
		o.type = CCT.COLLISION_BODY_SPHERE;
		if (pos)
			o.pos = new CCT.Vector3(pos.x, pos.y, pos.z);
		o.radius = radius;
		return o;
	},
	
	Capsule : function (pos, axis, radius, half_height) {
		var o = {};
		o.type= CCT.COLLISION_BODY_CAPSULE;
		if (pos)
			o.pos = new CCT.Vector3(pos.x, pos.y, pos.z);
		o.axis = new CCT.Vector3(axis.x, axis.y, axis.z);
		o.radius = radius;
		o.half_height = half_height;
		return o;
	},
	
	Plane : function (vertice, normal) {
		var o = {};
		o.type = CCT.COLLISION_BODY_PLANE;
		o.vertice = new CCT.Vector3(vertice.x, vertice.y, vertice.z);
		o.normal = new CCT.Vector3(normal.x, normal.y, normal.z);
		return o;
	},
	
	TrianglesPlane : function (normal, vertices, indices) {
		var o = {};
		o.type = CCT.COLLISION_BODY_TRIANGLES_PLANE;
		o.normal = new CCT.Vector3(normal.x, normal.y, normal.z);
		o.vertices = [].concat(vertices);
		o.indices = [].concat(indices);
		return o;
	},
	
	/**
	 * collision body1 cast body2
	 * @param one
	 * @param dir
	 * @param two
	 * @returns {object}
	 */
	 cast : function(one, dir, two) {
		if (one === two || CCT.mathVec3IsZero(dir))
			return null;
		else if (CCT.COLLISION_BODY_RAY === one.type) {
			switch (two.type) {
				case CCT.COLLISION_BODY_AABB:
					return CCT.mathRaycastAABB(one.pos, dir, two.pos, two.half);
				case CCT.COLLISION_BODY_SPHERE:
					return CCT.mathRaycastSphere(one.pos, dir, two.pos, two.radius);
				case CCT.COLLISION_BODY_CAPSULE:
					return CCT.mathRaycastCapsule(one.pos, dir, two.pos, two.axis, two.radius, two.half_height);
				case CCT.COLLISION_BODY_PLANE:
					return CCT.mathRaycastPlane(one.pos, dir, two.vertice, two.normal);
				case CCT.COLLISION_BODY_TRIANGLES_PLANE:
					return CCT.mathRaycastTrianglesPlane(one.pos, dir, two.normal, two.vertices, two.indices);
				default:
					return null;
			}
		}
		else if (CCT.COLLISION_BODY_AABB === one.type) {
			switch (two.type) {
				case CCT.COLLISION_BODY_AABB:
					return CCT.mathAABBcastAABB(one.pos, one.half, dir, two.pos, two.half);
				case CCT.COLLISION_BODY_SPHERE:
				{
					var result = CCT.mathSpherecastAABB(two.pos, two.radius, dir.clone().negate(), one.pos, one.half);
					if (result && result.hit_point)
						result.hit_point.addScaledVector(dir, result.distance);
					return result;
				}
				case CCT.COLLISION_BODY_CAPSULE:
				{
					var result = CCT.mathCapsulecastAABB(two.pos, two.axis, two.radius, two.half_height, dir.clone().negate(), one.pos, one.half);
					if (result && result.hit_point)
						result.hit_point.addScaledVector(dir, result.distance);
					return result;
				}
				case CCT.COLLISION_BODY_PLANE:
					return CCT.mathAABBcastPlane(one.pos, one.half, dir, two.vertice, two.normal);
				default:
					return null;
			}
		}
		else if (CCT.COLLISION_BODY_SPHERE === one.type) {
			switch (two.type) {
				case CCT.COLLISION_BODY_AABB:
					return CCT.mathSpherecastAABB(one.pos, one.radius, dir, two.pos, two.half);
				case CCT.COLLISION_BODY_SPHERE:
					return CCT.mathSpherecastSphere(one.pos, one.radius, dir, two.pos, two.radius);
				case CCT.COLLISION_BODY_CAPSULE:
					return CCT.mathSpherecastCapsule(one.pos, one.radius, dir, two.pos, two.axis, two.radius, two.half_height);
				case CCT.COLLISION_BODY_PLANE:
					return CCT.mathSpherecastPlane(one.pos, one.radius, dir, two.vertice, two.normal);
				case CCT.COLLISION_BODY_TRIANGLES_PLANE:
					return CCT.mathSpherecastTrianglesPlane(one.pos, one.radius, dir, two.normal, two.vertices, two.indices);
				default:
					return null;
			}
		}
		else if (CCT.COLLISION_BODY_CAPSULE === one.type) {
			switch (two.type) {
				case CCT.COLLISION_BODY_AABB:
				{
					return CCT.mathCapsulecastAABB(one.pos, one.axis, one.radius, one.half_height, dir, two.pos, two.half);
				}
				case CCT.COLLISION_BODY_SPHERE:
				{
					var result = CCT.mathSpherecastCapsule(two.pos, two.radius, dir.clone().negate(), one.pos, one.axis, one.radius, one.half_height);
					if (result && result.hit_point)
						result.hit_point.addScaledVector(dir, result.distance);
					return result;
				}
				case CCT.COLLISION_BODY_CAPSULE:
					return CCT.mathCapsulecastCapsule(one.pos, one.axis, one.radius, one.half_height, dir, two.pos, two.axis, two.radius, two.half_height);
				case CCT.COLLISION_BODY_PLANE:
					return CCT.mathCapsulecastPlane(one.pos, one.axis, one.radius, one.half_height, dir, two.vertice, two.normal);
				case CCT.COLLISION_BODY_TRIANGLES_PLANE:
					return CCT.mathCapsulecastTrianglesPlane(one.pos, one.axis, one.radius, one.half_height, dir, two.normal, two.vertices, two.indices);
				default:
					return null;
			}
		}
		return null;
	 },

	/**
	 * a==b return =0, a<b return <0, a>b return >0
	 * @param a
	 * @param b
	 * @param epsilon
	 * @returns {number}
	 */
	fcmpf : function(a, b, epsilon) {
		var v = parseFloat(a) - parseFloat(b);
		/* a == b */
		if (v > -epsilon && v < epsilon)
			return 0;
		return v >= epsilon ? 1 : -1;
	},

	mathVec3IsZero : function(v) {
		return	CCT.fcmpf(v.x, 0.0, CCT.EPSILON) === 0 &&
				CCT.fcmpf(v.y, 0.0, CCT.EPSILON) === 0 &&
				CCT.fcmpf(v.z, 0.0, CCT.EPSILON) === 0;
	},

	mathVec3Equal : function(a, b) {
		return	CCT.fcmpf(a.x, b.x, CCT.EPSILON) === 0 &&
				CCT.fcmpf(a.y, b.y, CCT.EPSILON) === 0 &&
				CCT.fcmpf(a.z, b.z, CCT.EPSILON) === 0;
	},

	mathCoordinateSystemTransform : function (v, new_origin, new_axies) {
		var t;
		if (new_origin)/* if v is normal vector, this field must be NULL */
			t = v.clone().sub(new_origin);
		else
			t = v;
		return new CCT.Vector3(
			t.dot(new_axies[0]),
			t.dot(new_axies[1]),
			t.dot(new_axies[2])
		);
	},

	mathPlaneNormalByVertices3 : function (vertices) {
		var v0v1 = new CCT.Vector3().subVectors(vertices[1], vertices[0]);
		var v0v2 = new CCT.Vector3().subVectors(vertices[2], vertices[0]);
		var normal = new CCT.Vector3().crossVectors(v0v1, v0v2);
		if (CCT.mathVec3IsZero(normal))
			return normal;
		return normal.normalize();
	},
	
	mathQuadraticEquation : function (a, b, c, r) {
		if (CCT.fcmpf(a, 0.0, CCT.EPSILON) === 0)
			return 0;
		var delta = b * b - 4.0 * a * c;
		var cmp = CCT.fcmpf(delta, 0.0, CCT.EPSILON);
		if (cmp < 0)
			return 0;
		else if (0 === cmp) {
			r[0] = r[1] = -b / a * 0.5;
			return 1;
		}
		else {
			var sqrt_delta = Math.sqrt(delta);
			r[0] = (-b + sqrt_delta) / a * 0.5;
			r[1] = (-b - sqrt_delta) / a * 0.5;
			return 2;
		}
	},
	
	mathPointProjectionLine : function (p, ls_v, lsdir) {
		var vp = p.clone().sub(ls_v);
		var dot = vp.dot(lsdir);
		var pj_p = ls_v.clone().addScaledVector(lsdir, dot);
		return {
			p : pj_p,
			v : p.clone().sub(pj_p)
		};
	},

	mathLineClosestLine : function (lsv1, lsdir1, lsv2, lsdir2) {
		var v = lsv2.clone().sub(lsv1);
		var n = new CCT.Vector3().crossVectors(lsdir1, lsdir2);
		var dot;
		if (CCT.mathVec3IsZero(n)) {
			dot = v.dot(lsdir1);
			return {
				code : 0,
				min_d : Math.sqrt(v.lengthSq() - dot * dot)
			};
		}
		var nlensq_inv = 1.0 / n.lengthSq();
		var N = n.clone().normalize();
		dot = v.dot(N);
		if (dot <= -CCT.EPSILON)
			dot = -dot;
		return {
			code : 1,
			min_d : dot,
			dir_d : [
				new CCT.Vector3().crossVectors(v, lsdir2).dot(n) * nlensq_inv,
				new CCT.Vector3().crossVectors(v, lsdir1).dot(n) * nlensq_inv
			]
		};
	},

	mathPointProjectionPlane : function (p, plane_v, plane_n) {
		var pv = new CCT.Vector3().subVectors(plane_v, p);
		var distance = pv.dot(plane_n);
		var point = p.clone().addScaledVector(plane_n, distance);
		return { p: point, d: distance };
	},
	
	mathPlaneHasPoint : function (plane_v, plane_n, p) {
		var v = new CCT.Vector3().subVectors(plane_v, p);
		return CCT.fcmpf(plane_n.dot(v), 0.0, CCT.EPSILON) === 0;
	},

	mathSegmentHasPoint : function (ls, p) {
		var v1 = ls[0], v2 = ls[1];
		var pv1 = new CCT.Vector3().subVectors(v1, p);
		var pv2 = new CCT.Vector3().subVectors(v2, p);
		if (!CCT.mathVec3IsZero(new CCT.Vector3().crossVectors(pv1, pv2)))
			return 0;
		else if (CCT.mathVec3Equal(ls[0], p))
			return 1;
		else if (CCT.mathVec3Equal(ls[1], p))
			return 2;
		else {
			var min_x, max_x, min_y, max_y, min_z, max_z;
			v1.x < v2.x ? (min_x = v1.x, max_x = v2.x) : (min_x = v2.x, max_x = v1.x);
			if (CCT.fcmpf(p.x, min_x, CCT.EPSILON) < 0 || CCT.fcmpf(p.x, max_x, CCT.EPSILON) > 0)
				return 0;
			v1.y < v2.y ? (min_y = v1.y, max_y = v2.y) : (min_y = v2.y, max_y = v1.y);
			if (CCT.fcmpf(p.y, min_y, CCT.EPSILON) < 0 || CCT.fcmpf(p.y, max_y, CCT.EPSILON) > 0)
				return 0;
			v1.z < v2.z ? (min_z = v1.z, max_z = v2.z) : (min_z = v2.z, max_z = v1.z);
			if (CCT.fcmpf(p.z, min_z, CCT.EPSILON) < 0 || CCT.fcmpf(p.z, max_z, CCT.EPSILON) > 0)
				return 0;
			return 1;
		}
	},

	mathTriangleHasPoint : function (tri, p) {
		var ap = new CCT.Vector3().subVectors(p, tri[0]);
		var ab = new CCT.Vector3().subVectors(tri[1], tri[0]);
		var ac = new CCT.Vector3().subVectors(tri[2], tri[0]);
		var N = new CCT.Vector3().crossVectors(ab, ac);
		if (CCT.fcmpf(N.dot(ap), 0.0, CCT.EPSILON))
			return null;
		else {
			var dot_ac_ac = ac.dot(ac);
			var dot_ac_ab = ac.dot(ab);
			var dot_ac_ap = ac.dot(ap);
			var dot_ab_ab = ab.dot(ab);
			var dot_ab_ap = ab.dot(ap);
			var tmp = 1.0 / (dot_ac_ac * dot_ab_ab - dot_ac_ab * dot_ac_ab);
			var u = (dot_ab_ab * dot_ac_ap - dot_ac_ab * dot_ab_ap) * tmp;
			if (CCT.fcmpf(u, 0.0, CCT.EPSILON) < 0 || CCT.fcmpf(u, 1.0, CCT.EPSILON) > 0)
				return null;
			var v = (dot_ac_ac * dot_ab_ap - dot_ac_ab * dot_ac_ap) * tmp;
			if (CCT.fcmpf(v, 0.0, CCT.EPSILON) < 0 || CCT.fcmpf(v + u, 1.0, CCT.EPSILON) > 0)
				return null;
			return { u: u, v: v};
		}
	},

	mathSphereHasPoint : function (o, radius, p) {
		var op = new CCT.Vector3().subVectors(p, o);
		var cmp = CCT.fcmpf(op.lengthSq(), radius * radius, CCT.EPSILON);
		if (cmp > 0)
			return 0;
		if (0 === cmp)
			return 1;
		else
			return 2;
	},
	
	mathCapsuleHasPoint : function (o, axis, radius, half_height, p) {
		var cp = o.clone().addScaledVector(axis, -half_height);
		var v = p.clone().sub(cp);
		var dot = axis.dot(v);
		if (CCT.fcmpf(dot, 0.0, CCT.EPSILON) < 0) {
			return CCT.mathSphereHasPoint(cp, radius, p);
		}
		else if (CCT.fcmpf(dot, half_height + half_height, CCT.EPSILON) > 0) {
			cp.addScaledVector(axis, half_height + half_height);
			return CCT.mathSphereHasPoint(cp, radius, p);
		}
		else {
			var cmp = CCT.fcmpf(v.lengthSq() - dot * dot, radius * radius, CCT.EPSILON);
			if (cmp > 0)
				return 0;
			return cmp < 0 ? 2 : 1;
		}
	},
	
	mathAABBHasPoint : function (o, half, p) {
		return p.x >= o.x - half.x && p.x <= o.x + half.x &&
			p.y >= o.y - half.y && p.y <= o.y + half.y &&
			p.z >= o.z - half.z && p.z <= o.z + half.z;
	},
	
	mathLineIntersectLine : function (ls1v, ls1dir, ls2v, ls2dir) {
		var dot, N = new CCT.Vector3().crossVectors(ls1dir, ls2dir);
		var v = ls1v.clone().sub(ls2v);
		if (CCT.mathVec3IsZero(N)) {
			dot = v.dot(ls2dir);
			return CCT.fcmpf(dot * dot, v.lengthSq(), CCT.EPSILON) ? null : { code : 2 };
		}
		else if (CCT.mathVec3IsZero(v)) {
			return {
				code : 1,
				dir_d : [ 0.0, 0.0 ]
			};
		}
		dot = v.dot(N);
		if (CCT.fcmpf(dot, 0.0, CCT.EPSILON))
			return null;
		dot = v.dot(ls2dir);
		v = ls2v.clone().addScaledVector(ls2dir, dot).sub(ls1v);
		var dir_d = new Array(2);
		if (CCT.mathVec3IsZero(v)) {
			dir_d[0] = 0.0;
		}
		else {
			dir_d[0] = v.length();
			dir_d[0] /= v.normalize().dot(ls1dir);
		}
		v = ls2v.clone().sub(ls1v);
		dot = v.dot(ls1dir);
		v = ls1v.clone().addScaledVector(ls1dir, dot).sub(ls2v);
		if (CCT.mathVec3IsZero(v)) {
			dir_d[1] = 0.0;
		}
		else {
			dir_d[1] = v.length();
			dir_d[1] /= v.normalize().dot(ls2dir);
		}
		return {
			code : 1,
			dir_d : dir_d
		};
	},
	
	overlapSegmentIntersectSegment : function (ls1, ls2) {
		var res = CCT.mathSegmentHasPoint(ls1, ls2[0]);
		if (3 === res)
			return { code : 2 };
		else if (res) {
			if (CCT.mathSegmentHasPoint(ls1, ls2[1]))
				return { code : 2 };
			if (CCT.mathSegmentHasPoint(ls2, ls1[res === 1 ? 1 : 0]))
				return { code : 2 };
			return {
				code : 1,
				p : ls1[res - 1].clone()
			};
		}
		res = CCT.mathSegmentHasPoint(ls1, ls2[1]);
		if (3 === res)
			return { code : 2 };
		else if (res) {
			if (CCT.mathSegmentHasPoint(ls2, ls1[res === 1 ? 1 : 0]))
				return { code : 2 };
			return {
				code : 1,
				p : ls1[res - 1].clone()
			};
		}

		if (CCT.mathSegmentHasPoint(ls2, ls1[0]) === 3 ||
			CCT.mathSegmentHasPoint(ls2, ls1[1]) === 3)
		{
			return { code : 2 };
		}
		return null;
	},
	
	mathSegmentIntersectSegment : function (ls1, ls2) {
		var lsdir1 = ls1[1].clone().sub(ls1[0]);
		var lslen1 = lsdir1.length();
		lsdir1.normalize();
		var lsdir2 = ls2[1].clone().sub(ls2[0]);
		var lslen2 = lsdir2.length();
		lsdir2.normalize();
		var res = CCT.mathLineIntersectLine(ls1[0], lsdir1, ls2[0], lsdir2);
		if (!res)
			return null;
		else if (1 === res.code) {
			if (CCT.fcmpf(res.dir_d[0], 0.0, CCT.EPSILON) < 0 || CCT.fcmpf(res.dir_d[1], 0.0, CCT.EPSILON) < 0)
				return null;
			if (CCT.fcmpf(res.dir_d[0], lslen1, CCT.EPSILON) > 0 || CCT.fcmpf(res.dir_d[1], lslen2, CCT.EPSILON) > 0)
				return null;
			return {
				code : 1,
				p : ls1[0].clone().addScaledVector(lsdir1, res.dir_d[0])
			};
		}
		else
			return CCT.overlapSegmentIntersectSegment(ls1, ls2);
	},
	
	mathLineIntersectPlane : function (ls_v, lsdir, plane_v, plane_n) {
		var cos_theta = lsdir.dot(plane_n);
		var res = CCT.mathPointProjectionPlane(ls_v, plane_v, plane_n);
		if (CCT.fcmpf(cos_theta, 0.0, CCT.EPSILON)) {
			return {
				code : 1,
				dir_d : res.d / cos_theta
			};
		}
		else {
			return CCT.fcmpf(res.d, 0.0, CCT.EPSILON) ? null : { code : 2 };
		}
	},
	
	mathSegmentIntersectPlane : function (ls, plane_v, plane_n) {
		var pj = [
			CCT.mathPointProjectionPlane(ls[0], plane_v, plane_n),
			CCT.mathPointProjectionPlane(ls[1], plane_v, plane_n)
		];
		var cmp = [
			CCT.fcmpf(pj[0].d, 0.0, CCT.EPSILON),
			CCT.fcmpf(pj[1].d, 0.0, CCT.EPSILON)
		];
		if (0 === cmp[0] && 0 === cmp[1])
			return { code : 2 };
		else if (cmp[0] * cmp[1] > 0)
			return null;
		else if (0 === cmp[0]) {
			return {
				code : 1,
				p : ls[0].clone()
			};
		}
		else if (0 === cmp[1]) {
			return {
				code : 1,
				p : ls[1].clone()
			};
		}
		else {
			var lsdir = ls[1].clone().sub(ls[0]).normalize();
			var dot = lsdir.dot(plane_n);
			return {
				code : 1,
				p : ls[0].clone().addScaledVector(lsdir, pj[0].d / dot)
			};
		}
	},
	
	mathSphereIntersectLine : function (o, radius, ls_vertice, lsdir) {
		var vo = o.clone().sub(ls_vertice);
		var dot = vo.dot(lsdir);
		var lp = ls_vertice.clone().addScaledVector(lsdir, dot);
		var lpo = o.clone().sub(lp);
		var lpolensq = lpo.lengthSq();
		var radiussq = radius * radius;
		var cmp = CCT.fcmpf(lpolensq, radiussq, CCT.EPSILON);
		if (cmp > 0)
			return null;
		else if (0 === cmp) {
			return {
				code : 1,
				dir_d : [ dot, dot ]
			};
		}
		else {
			var d = Math.sqrt(radiussq - lpolensq);
			return {
				code : 2,
				dir_d : [ dot + d, dot - d ]
			};
		}
	},
	
	mathSphereIntersectSegment : function (o, radius, ls) {
		var c = [
			CCT.mathSphereHasPoint(o, radius, ls[0]),
			CCT.mathSphereHasPoint(o, radius, ls[1])
		];
		if (c[0] + c[1] >= 2)
			return { code : 2 };
		else {
			var lsdir = ls[1].clone().sub(ls[0]).normalize();
			var pj = CCT.mathPointProjectionLine(o, ls[0], lsdir);
			if (!CCT.mathSegmentHasPoint(ls, pj.p)) {
				if (0 === c[0] + c[0])
					return null;
				var res = { code : c[0] + c[1] };
				if (c[0])
					res.p = ls[0].clone();
				else if (c[1])
					res.p = ls[1].clone();
				return res;
			}
			c[0] = CCT.fcmpf(pj.v.lengthSq(), radius * radius, CCT.EPSILON);
			if (c[0] < 0)
				return { code : 2 };
			if (0 === c[0]) {
				return {
					code : 1,
					p : pj.p
				};
			}
			return null;
		}
	},
	
	mathSphereIntersectPlane : function (o, radius, plane_v, plane_n) {
		var pj = CCT.mathPointProjectionPlane(o, plane_v, plane_n);
		var ppo = o.clone().sub(pj.p);
		var ppolensq = ppo.lengthSq();
		var cmp = CCT.fcmpf(ppolensq, radius * radius, CCT.EPSILON);
		if (cmp > 0)
			return null;
		var res = { p : pj.p };
		if (0 === cmp) {
			res.code = 1;
			res.r = 0.0;
		}
		else {
			res.code = 2;
			res.r = Math.sqrt(radius * radius - ppolensq);
		}
		return res;
	},
	
	mathSphereIntersectTrianglesPlane : function (o, radius, plane_n, vertices, indices) {
		var res = CCT.mathSphereIntersectPlane(o, radius, vertices[indices[0]], plane_n);
		if (!res)
			return null;
		var i = 0;
		while (i < indices.length) {
			var tri = [
				vertices[indices[i++]],
				vertices[indices[i++]],
				vertices[indices[i++]]
			];
			if (CCT.mathTriangleHasPoint(tri, res.p))
				return { code : 1 };
		}
		if (res.code === 2) {
			for (i = 0; i < indices.length; i += 3) {
				for (var j = 0; j < 3; ++j) {
					var edge = [
						vertices[indices[j % 3 + i]],
						vertices[indices[(j + 1) % 3 + i]]
					];
					if (CCT.mathSphereIntersectSegment(res.p, res.r, edge))
						return { code : 1 };
				}
			}
		}
		return null;
	},
	
	mathSphereIntersectSphere : function (o1, r1, o2, r2) {
		var radius_sum = r1 + r2;
		var o1o2 = o2.clone().sub(o1);
		var cmp = CCT.fcmpf(o1o2.lengthSq(), radius_sum * radius_sum, CCT.EPSILON);
		if (cmp > 0)
			return null;
		if (cmp < 0)
			return { code : 2 };
		else {
			return {
				code : 1,
				p : o1.clone().addScaledVector(o1o2.normalize(), r1)
			};
		}
	},
	
	mathSphereIntersectCapsule : function (sp_o, sp_radius, cp_o, cp_axis, cp_radius, cp_half_height) {
		var cp = cp_o.clone().addScaledVector(cp_axis, -cp_half_height);
		var v = sp_o.clone().sub(cp);
		var dot = cp_axis.dot(v);
		if (CCT.fcmpf(dot, 0.0, CCT.EPSILON) < 0) {
			return CCT.mathSphereIntersectSphere(sp_o, sp_radius, cp, cp_radius);
		}
		else if (CCT.fcmpf(dot, cp_half_height + cp_half_height, CCT.EPSILON) > 0) {
			cp.addScaledVector(cp_axis, cp_half_height + cp_half_height);
			return CCT.mathSphereIntersectSphere(sp_o, sp_radius, cp, cp_radius);
		}
		else {
			var radius_sum = sp_radius + cp_radius;
			var cmp = CCT.fcmpf(v.lengthSq() - dot * dot, radius_sum * radius_sum, CCT.EPSILON);
			if (cmp > 0)
				return null;
			if (cmp < 0)
				return { code : 2 };
			else {
				v = cp.clone().addScaledVector(cp_axis, dot).sub(sp_o).normalize();
				return {
					code : 1,
					p : sp_o.clone().addScaledVector(v, sp_radius)
				};
			}
		}
	},
	
	Box_Edge_Indices : [
		0, 1,	1, 2,	2, 3,	3, 0,
		7, 6,	6, 5,	5, 4,	4, 7,
		1, 5,	6, 2,
		3, 7,	4, 0
	],

	Box_Triangle_Vertices_Indices : [
		0, 1, 2,	2, 3, 0,
		7, 6, 5,	5, 4, 7,
		1, 5, 6,	6, 2, 1,
		3, 7, 4,	4, 0, 3,
		3, 7, 6,	6, 2, 3,
		0, 4, 5,	5, 1, 0
	],
	
	AABB_Plane_Normal : [
		[ 0.0, 0.0, 1.0 ], [ 0.0, 0.0, -1.0 ],
		[ 1.0, 0.0, 0.0 ], [ -1.0, 0.0, 0.0 ],
		[ 0.0, 1.0, 0.0 ], [ 0.0, -1.0, 0.0 ]
	],
	
	AABBVertices : function(o, half) {
		return [
			new CCT.Vector3(o.x - half.x, o.y - half.y, o.z - half.z),
			new CCT.Vector3(o.x + half.x, o.y - half.y, o.z - half.z),
			new CCT.Vector3(o.x + half.x, o.y + half.y, o.z - half.z),
			new CCT.Vector3(o.x - half.x, o.y + half.y, o.z - half.z),
			new CCT.Vector3(o.x - half.x, o.y - half.y, o.z + half.z),
			new CCT.Vector3(o.x + half.x, o.y - half.y, o.z + half.z),
			new CCT.Vector3(o.x + half.x, o.y + half.y, o.z + half.z),
			new CCT.Vector3(o.x - half.x, o.y + half.y, o.z + half.z)
		];
	},
	
	mathAABBIntersectAABB : function (o1, half1, o2, half2) {
		return !(o2.x - o1.x > half1.x + half2.x || o1.x - o2.x > half1.x + half2.x ||
				o2.y - o1.y > half1.y + half2.y || o1.y - o2.y > half1.y + half2.y ||
				o2.z - o1.z > half1.z + half2.z || o1.z - o2.z > half1.z + half2.z);
	},
	
	mathAABBIntersectPlane : function (o, half, plane_v, plane_n) {
		var prev_d;
		var vertices = CCT.AABBVertices(o, half, vertices);
		for (var i = 0; i < 8; ++i) {
			var pj = CCT.mathPointProjectionPlane(vertices[i], plane_v, plane_n);
			if (i && prev_d * pj.d <= -CCT.EPSILON)
				return null;
			prev_d = pj.d;
		}
		return { code : 1 };
	},
	
	mathAABBIntersectSphere : function (aabb_o, aabb_half, sp_o, sp_radius) {
		if (CCT.mathAABBHasPoint(aabb_o, aabb_half, sp_o) || CCT.mathSphereHasPoint(sp_o, sp_radius, aabb_o))
			return { code : 1 };
		else {
			var v = CCT.AABBVertices(aabb_o, aabb_half);
			for (var i = 0, j = 0; i < CCT.Box_Triangle_Vertices_Indices.length; i += 6, ++j) {
				var plane_n = new CCT.Vector3().fromArray(CCT.AABB_Plane_Normal[j]);
				if (CCT.mathSphereIntersectTrianglesPlane(sp_o, sp_radius,
					plane_n, v, CCT.Box_Triangle_Vertices_Indices.slice(i, i + 6)))
				{
					return {code: 1};
				}
			}
			return null;
		}
	},
	
	mathAABBIntersectCapsule : function (aabb_o, aabb_half, cp_o, cp_axis, cp_radius, cp_half_height) {
		if (CCT.mathAABBHasPoint(aabb_o, aabb_half, cp_o) ||
			CCT.mathCapsuleHasPoint(cp_o, cp_axis, cp_radius, cp_half_height, aabb_o))
		{
			return { code: 1 };
		}
		else {
			var v = CCT.AABBVertices(aabb_o, aabb_half);
			for (var i = 0, j = 0; i < CCT.Box_Triangle_Vertices_Indices.length; i += 6, ++j) {
				var plane_n = new CCT.Vector3().fromArray(CCT.AABB_Plane_Normal[j]);
				if (CCT.mathCapsuleIntersectTrianglesPlane(cp_o, cp_axis, cp_radius, cp_half_height,
					plane_n, v, CCT.Box_Triangle_Vertices_Indices.slice(i, i + 6)))
				{
					return {code: 1};
				}
			}
			return null;
		}
	},
	
	mathLineIntersectCylinderInfinite : function (ls_v, lsdir, cp, axis, radius) {
		var radius_sq = radius * radius;
		var new_axies = new Array(3);
		new_axies[2] = axis.clone();
		new_axies[1] = new CCT.Vector3(0.0, -new_axies[2].z, new_axies[2].y);
		if (CCT.mathVec3IsZero(new_axies[1])) {
			new_axies[1] = new CCT.Vector3(-new_axies[2].z, 0.0, new_axies[2].x);
		}
		new_axies[1].normalize();
		new_axies[0] = new CCT.Vector3().crossVectors(new_axies[1], new_axies[2]);
		var new_o = CCT.mathCoordinateSystemTransform(ls_v, cp, new_axies);
		var new_dir = CCT.mathCoordinateSystemTransform(lsdir, null, new_axies);
		var A = new_dir.x * new_dir.x + new_dir.y * new_dir.y;
		var B = 2.0 * (new_o.x * new_dir.x + new_o.y * new_dir.y);
		var C = new_o.x * new_o.x + new_o.y * new_o.y - radius_sq;
		var r = new Array(2);
		var rcnt = CCT.mathQuadraticEquation(A, B, C, r);
		if (0 === rcnt) {
			var v = ls_v.clone().sub(cp);
			var dot = v.dot(axis);
			return CCT.fcmpf(v.lengthSq() - dot * dot, radius_sq, CCT.EPSILON) > 0 ? null : { code : -1 };
		}
		return {
			code : rcnt,
			dir_d : r
		};
	},
	
	mathCapsuleIntersectPlane : function (cp_o, cp_axis, cp_radius, cp_half_height, plane_v, plane_n) {
		var sphere_o = new Array(2), pj = new Array(2), cmp = new Array(2);
		for (var i = 0; i < 2; ++i) {
			sphere_o[i] = cp_o.clone().addScaledVector(cp_axis, i ? cp_half_height : -cp_half_height);
			pj[i] = CCT.mathPointProjectionPlane(sphere_o[i], plane_v, plane_n);
			cmp[i] = CCT.fcmpf(pj[i].d * pj[i].d, cp_radius * cp_radius, CCT.EPSILON);
			if (cmp[i] < 0)
				return { code : 2 };
		}
		if (0 === cmp[0] && 0 === cmp[1])
			return { code : 2 };
		else if (0 === cmp[0]) {
			return {
				code : 1,
				p : sphere_o[0].clone().addScaledVector(plane_n, pj[0].d)
			};
		}
		else if (0 === cmp[1]) {
			return {
				code : 1,
				p : sphere_o[1].clone().addScaledVector(plane_n, pj[1].d)
			};
		}
		cmp[0] = CCT.fcmpf(pj[0].d, 0.0, CCT.EPSILON);
		cmp[1] = CCT.fcmpf(pj[1].d, 0.0, CCT.EPSILON);
		return cmp[0] * cmp[1] > 0 ? null : { code : 2 };
	},
	
	mathCapsuleIntersectTrianglesPlane : function (cp_o, cp_axis, cp_radius, cp_half_height, plane_n, vertices, indices) {
		var i, res = CCT.mathCapsuleIntersectPlane(cp_o, cp_axis, cp_radius, cp_half_height, vertices[indices[0]], plane_n);
		if (!res)
			return null;
		else if (1 === res) {
			i = 0;
			while (i < indices.length) {
				var tri = [
					vertices[indices[i++]],
					vertices[indices[i++]],
					vertices[indices[i++]]
				];
				if (CCT.mathTriangleHasPoint(tri, res.p))
					return { code : 1 };
			}
			return null;
		}
		else {
			for (i = 0; i < indices.length; i += 3) {
				for (var j = 0; j < 3; ++j) {
					var edge = [
						vertices[indices[j % 3 + i]],
						vertices[indices[(j + 1) % 3 + i]]
					];
					if (CCT.mathSegmentIntersectCapsule(edge, cp_o, cp_axis, cp_radius, cp_half_height))
						return { code : 1 };
				}
			}
			var cos_theta = cp_axis.dot(plane_n);
			var center;
			if (CCT.fcmpf(cos_theta, 0.0, CCT.EPSILON)) {
				var pj = CCT.mathPointProjectionPlane(cp_o, vertices[indices[0]], plane_n);
				center = cp_o.clone().addScaledVector(cp_axis, pj.d);
			}
			else {
				center = CCT.mathPointProjectionPlane(cp_o, vertices[indices[0]], plane_n).p;
			}
			for (i = 0; i < indices.length; ) {
				var tri = [
					vertices[indices[i++]],
					vertices[indices[i++]],
					vertices[indices[i++]]
				];
				if (CCT.mathTriangleHasPoint(tri, center))
					return { code : 1 };
			}
			return null;
		}
	},
	
	mathCapsuleIntersectCapsule : function (cp1_o, cp1_axis, cp1_radius, cp1_half_height, cp2_o, cp2_axis, cp2_radius, cp2_half_height) {
		var radius_sum = cp1_radius + cp2_radius;
		var res = CCT.mathLineClosestLine(cp1_o, cp1_axis, cp2_o, cp2_axis);
		if (CCT.fcmpf(res.min_d, radius_sum, CCT.EPSILON) > 0)
			return null;
		else {
			var i, sphere_o;
			for (i = 0; i < 2; ++i) {
				sphere_o = cp1_o.clone().addScaledVector(cp1_axis, i ? cp1_half_height : -cp1_half_height);
				if (CCT.mathSphereIntersectCapsule(sphere_o, cp1_radius, cp2_o, cp2_axis, cp2_radius, cp2_half_height))
					return { code : 1 };
			}
			for (i = 0; i < 2; ++i) {
				sphere_o = cp2_o.clone().addScaledVector(cp2_axis, i ? cp2_half_height : -cp2_half_height);
				if (CCT.mathSphereIntersectCapsule(sphere_o, cp2_radius, cp1_o, cp1_axis, cp1_radius, cp1_half_height))
					return { code : 1 };
			}
			if (0 === res.code)
				return null;
			if (CCT.fcmpf(res.dir_d[0] >= CCT.EPSILON ? res.dir_d[0] : -res.dir_d[0], cp1_half_height, CCT.EPSILON) < 0 &&
				CCT.fcmpf(res.dir_d[1] >= CCT.EPSILON ? res.dir_d[1] : -res.dir_d[1], cp2_half_height, CCT.EPSILON) < 0)
			{
				return { code : 1 };
			}
			return null;
		}
	},
	
	mathLineIntersectCapsule : function (ls_v, lsdir, o, axis, radius, half_height) {
		var res = CCT.mathLineIntersectCylinderInfinite(ls_v, lsdir, o, axis, radius);
		if (!res)
			return null;
		if (1 === res.code) {
			var p = ls_v.clone().addScaledVector(lsdir, res.dir_d[0]);
			var code = CCT.mathCapsuleHasPoint(o, axis, radius, half_height, p);
			if (code) {
				return {
					code : code,
					dir_d : res.dir_d
				};
			}
			return null;
		}
		else {
			var i, j = -1, cnt;
			if (2 === res.code) {
				cnt = 0;
				for (i = 0; i < 2; ++i) {
					var p = ls_v.clone().addScaledVector(lsdir, res.dir_d[i]);
					if (!CCT.mathCapsuleHasPoint(o, axis, radius, half_height, p))
						continue;
					++cnt;
					j = i;
				}
				if (2 === cnt)
					return { code : 2, dir_d : res.dir_d };
			}
			cnt = 0;
			var sphere_o, obj, d = new Array(5);
			sphere_o = o.clone().addScaledVector(axis, half_height);
			obj = CCT.mathSphereIntersectLine(sphere_o, radius, ls_v, lsdir);
			if (obj) {
				d[cnt++] = obj.dir_d[0];
				d[cnt++] = obj.dir_d[1];
			}
			sphere_o = o.clone().addScaledVector(axis, -half_height);
			obj = CCT.mathSphereIntersectLine(sphere_o, radius, ls_v, lsdir);
			if (obj) {
				d[cnt++] = obj.dir_d[0];
				d[cnt++] = obj.dir_d[1];
			}
			if (0 === cnt)
				return null;
			var min_d = null, max_d = null;
			if (j >= 0)
				d[cnt++] = res.dir_d[j];
			for (i = 0; i < cnt; ++i) {
				if (null === min_d || min_d > d[i])
					min_d = d[i];
				if (null === max_d || max_d < d[i])
					max_d = d[i];
			}
			return {
				code : 2,
				dir_d : [ min_d, max_d ]
			};
		}
	},
	
	mathSegmentIntersectCapsule : function (ls, o, axis, radius, half_height) {
		var res = new Array(2);
		res[0] = CCT.mathCapsuleHasPoint(o, axis, radius, half_height, ls[0]);
		if (2 === res[0])
			return { code : 2 };
		res[1] = CCT.mathCapsuleHasPoint(o, axis, radius, half_height, ls[1]);
		if (2 === res[1])
			return { code : 2 };
		if (res[0] + res[1] >= 2)
			return { code : 2 };
		else {
			var lsdir = ls[1].clone().sub(ls[0]).normalize();
			res = CCT.mathLineIntersectCapsule(ls[0], lsdir, o, axis, radius, half_height);
			if (!res)
				return null;
			var new_ls = [
				ls[0].clone().addScaledVector(lsdir, res.dir_d[0]),
				ls[0].clone().addScaledVector(lsdir, res.dir_d[1])
			];
			return CCT.overlapSegmentIntersectSegment(ls, new_ls);
		}
	},

	mathRaycastSegment : function (o, dir, ls) {
		var result;
		var lsdir = ls[1].clone().sub(ls[0]);
		var lslen = lsdir.length();
		lsdir.normalize();
		var res = CCT.mathLineIntersectLine(o, dir, ls[0], lsdir);
		if (!res)
			return null;
		else if (1 === res.code) {
			if (CCT.fcmpf(res.dir_d[0], 0.0, CCT.EPSILON) < 0)
				return null;
			if (CCT.fcmpf(res.dir_d[1], 0.0, CCT.EPSILON) < 0 || CCT.fcmpf(res.dir_d[1], lslen, CCT.EPSILON) > 0)
				return null;
			result = {};
			result.distance = res.dir_d[0];
			result.hit_point = o.clone().addScaledVector(dir, res.dir_d[0]);
			result.hit_normal = CCT.mathPointProjectionLine(o, ls[0], lsdir).v;
			return result;
		}
		else {
			result = {};
			result.hit_normal = ls[0].clone().sub(o);
			var d = [
				ls[0].clone().sub(o).dot(dir),
				ls[1].clone().sub(o).dot(dir)
			];
			var cmp = [
				CCT.fcmpf(d[0], 0.0, CCT.EPSILON),
				CCT.fcmpf(d[1], 0.0, CCT.EPSILON)
			];
			if (cmp[0] < 0 && cmp[0] < 0)
				return null;
			if (cmp[0] > 0 && cmp[0] > 0) {
				if (d[0] < d[1]) {
					result.distance = d[0];
					result.hit_point = ls[0].clone();
				}
				else {
					result.distance = d[1];
					result.hit_point = ls[1].clone();
				}
			}
			else {
				result.hit_point = o.clone();
			}
			return result;
		}
	},
	/**
	 * ray cast to plane
	 * @param o
	 * @param dir
	 * @param plane_v
	 * @param plane_n
	 * @returns {*}
	 */
	mathRaycastPlane : function (o, dir, plane_v, plane_n) {
		var pp = CCT.mathPointProjectionPlane(o, plane_v, plane_n);
		if (CCT.fcmpf(pp.d, 0.0, CCT.EPSILON) === 0) {
			return {
				distance : 0.0,
				hit_point : o.clone()
			};
		}
		var cos_theta = dir.dot(N);
		if (CCT.fcmpf(cos_theta, 0.0, CCT.EPSILON) === 0)
			return null;
		var distance = pp.d / cos_theta;
		if (CCT.fcmpf(distance, 0.0, CCT.EPSILON) < 0)
			return null;
		return {
			distance : distance,
			hit_point : o.clone().addScaledVector(dir, distance),
			hit_normal : plane_n.clone()
		};
	},
	/**
	 * ray cast to triangle
	 * @param o
	 * @param dir
	 * @param tri
	 * @returns {*}
	 */
	mathRaycastTriangle : function (o, dir, tri) {
		var result;
		if (result = CCT.mathRaycastPlane(o, dir, tri[0], CCT.mathPlaneNormalByVertices3(tri))) {
			if (CCT.mathTriangleHasPoint(tri, result.hit_point))
				return result;
			else if (CCT.fcmpf(result.distance, 0.0, CCT.EPSILON) === 0) {
				result = null;
				var results = new Array(3);
				for (var i = 0; i < 3; ++i) {
					var edge = [
						tri[i % 3],
						tri[(i + 1) % 3]
					];
					if (!(results[i] = CCT.mathRaycastSegment(o, dir, edge)))
						continue;
					if (!result || result.distance > results[i].distance) {
						result = results[i];
					}
				}
				return result;
			}
		}
		return null;
	},
	/**
	 * ray cast to triangles plane
	 * @param o
	 * @param dir
	 * @param plane_n
	 * @param vertices
	 * @param indices
	 * @returns {*}
	 */
	mathRaycastTrianglesPlane : function (o, dir, plane_n, vertices, indices) {
		var result = CCT.mathRaycastPlane(o, dir, vertices[indices[0]], plane_n);
		if (result) {
			var i = 0;
			while (i < indices.length) {
				var tri = [
					vertices[indices[i++]],
					vertices[indices[i++]],
					vertices[indices[i++]]
				];
				if (CCT.mathTriangleHasPoint(tri, result.hit_point))
					return result;
			}
		}
		return null;
	},
	/**
	 * ray cast to triangles plane
	 * @param o
	 * @param dir
	 * @param aabb_o
	 * @param aabb_half
	 * @returns {*}
	 */
	mathRaycastAABB : function (o, dir, aabb_o, aabb_half) {
		if (CCT.mathAABBHasPoint(aabb_o, aabb_half, o)) {
			return {
				distance : 0.0,
				hit_point : o.clone()
			};
		}
		else {
			var result = null;
			var v = CCT.AABBVertices(aabb_o, aabb_half);
			for (var i = 0, j = 0; i < CCT.Box_Triangle_Vertices_Indices.length; i += 6, ++j) {
				var plane_n = new CCT.Vector3().fromArray(CCT.AABB_Plane_Normal[j]);
				var result_temp = CCT.mathRaycastTrianglesPlane(o, dir,
					plane_n, v, CCT.Box_Triangle_Vertices_Indices.slice(i, i + 6));
				if (!result_temp)
					continue;
				if (!result || result.distance > result_temp.distance)
				{
					result = result_temp;
				}
			}
			return result;
		}
	},
	/**
	 * ray cast to sphere
	 * @param o
	 * @param dir
	 * @param center
	 * @param radius
	 * @returns {*}
	 */
	mathRaycastSphere : function (o, dir, center, radius) {
		var radius2 = radius * radius;
		var oc = new CCT.Vector3().subVectors(center, o);
		var oc2 = oc.lengthSq();
		var dir_d = dir.dot(oc);
		if (CCT.fcmpf(oc2, radius2, CCT.EPSILON) <= 0) {
			return {
				distance : 0.0,
				hit_point : o.clone()
			};
		}
		else if (CCT.fcmpf(dir_d, 0.0, CCT.EPSILON) <= 0)
			return null;

		var dr2 = oc2 - dir_d * dir_d;
		if (CCT.fcmpf(dr2, radius2, CCT.EPSILON) > 0)
			return null;

		var d = Math.sqrt(radius2 - dr2);
		var distance = dir_d - d;
		var result = {
			distance : dir_d - d,
			hit_point : o.clone().addScaledVector(dir, distance)
		};
		result.hit_normal = result.hit_point.clone().sub(center);
		return result;
	},
	mathRaycastCapsule : function (o, dir, cp_o, cp_axis, cp_radius, cp_half_height) {
		if (CCT.mathCapsuleHasPoint(cp_o, cp_axis, cp_radius, cp_half_height, o)) {
			return {
				distance : 0.0,
				hit_point : o.clone()
			};
		}
		else {
			var res;
			if (res = CCT.mathLineIntersectCapsule(o, dir, cp_o, cp_axis, cp_radius, cp_half_height)) {
				var cmp = [
					CCT.fcmpf(res.dir_d[0], 0.0, CCT.EPSILON),
					CCT.fcmpf(res.dir_d[1], 0.0, CCT.EPSILON)
				];
				if (cmp[0] > 0 && cmp[1] > 0) {
					var result = {};
					result.distance = res.dir_d[0] < res.dir_d[1] ? res.dir_d[0] : res.dir_d[1];
					result.hit_point = o.clone().addScaledVector(dir, result.distance);

					var sphere_o = cp_o.clone().addScaledVector(cp_axis, -cp_half_height);
					var v = result.hit_point.clone().sub(sphere_o);
					var dot = v.dot(cp_axis);
					if (CCT.fcmpf(dot, 0.0, CCT.EPSILON) < 0) {
						result.hit_normal = v;
					}
					else if (CCT.fcmpf(dot, cp_half_height, CCT.EPSILON) > 0) {
						sphere_o = cp_o.clone().addScaledVector(cp_axis, cp_half_height);
						v = result.hit_point.clone().sub(sphere_o);
						result.hit_normal = v;
					}
					else {
						result.hit_normal = CCT.mathPointProjectionLine(result.hit_point, cp_o, cp_axis).v;
					}
					return result;
				}
			}
			return null;
		}	
	},
	/**
	 * line segment cast to plane
	 * @param ls
	 * @param dir
	 * @param vertice
	 * @param normal
	 * @returns {*}
	 */
	mathSegmentcastPlane : function (ls, dir, vertice, normal) {
		var res = CCT.mathSegmentIntersectPlane(ls, vertice, normal);
		if (res && 2 === res.code) {
			return { distance : 0.0 };
		}
		else if (res && 1 === res.code) {
			return { distance : 0.0, hit_point : res.p };
		}
		else {
			var cos_theta = normal.dot(dir);
			if (CCT.fcmpf(cos_theta, 0.0, CCT.EPSILON) === 0)
				return null;
			var d = [
				CCT.mathPointProjectionPlane(ls[0], vertice, normal).d,
				CCT.mathPointProjectionPlane(ls[1], vertice, normal).d
			];
			var min_d, p, result = {};
			if (CCT.fcmpf(d[0], d[1], CCT.EPSILON) === 0) {
				min_d = d[0];
			}
			else {
				result.hit_point = {};
				if (CCT.fcmpf(d[0], 0.0, CCT.EPSILON) > 0) {
					if (d[0] < d[1]) {
						min_d = d[0];
						p = ls[0];
					}
					else {
						min_d = d[1];
						p = ls[1];
					}
				}
				else {
					if (d[0] < d[1]) {
						min_d = d[1];
						p = ls[1];
					}
					else {
						min_d = d[0];
						p = ls[0];
					}
				}
			}
			min_d /= cos_theta;
			if (CCT.fcmpf(min_d, 0.0, CCT.EPSILON) < 0)
				return null;
			result.distance = min_d;
			if (result.hit_point)
				result.hit_point = p.clone().addScaledVector(dir, result.distance);
			result.hit_normal = normal.clone();
			return result;
		}
	},
	/**
	 * line segment cast to line segment
	 * @param ls1
	 * @param dir
	 * @param ls2
	 * @returns {*}
	 */
	mathSegmentcastSegment : function (ls1, dir, ls2) {
		var res = CCT.mathSegmentIntersectSegment(ls1, ls2);
		if (res && 1 === res.code) {
			return {
				distance : 0.0,
				hit_point : res.p
			};
		}
		else if (res && 2 === res.code) {
			return { distance : 0.0 };
		}
		else {
			var lsdir1 = ls1[1].clone().sub(ls1[0]);
			var N = new CCT.Vector3().crossVectors(lsdir1, dir);
			if (CCT.mathVec3IsZero(N)) {
				var result0, result1;
				result0 = CCT.mathRaycastSegment(ls1[0], dir, ls2);
				if (!result0)
					return null;
				result1 = CCT.mathRaycastSegment(ls1[1], dir, ls2);
				if (!result1)
					return null;
				return result0.distance < result1.distance ? result0 : result1;
			}
			else {
				res = CCT.mathSegmentIntersectPlane(ls2, ls1[0], N.normalize());
				if (!res)
					return null;
				var neg_dir, result = null;
				if (1 === res.code) {
					var hit_point = res.p.clone();
					neg_dir = dir.clone().negate();
					result = CCT.mathRaycastSegment(hit_point, neg_dir, ls1);
					if (!result)
						return null;
					result.hit_point = hit_point;
					return result;
				}
				else if (2 === res.code) {
					var lsdir2 = ls2[1].clone().sub(ls2[0]);
					var v = new CCT.Vector3().crossVectors(lsdir1, lsdir2);
					var is_parallel = CCT.mathVec3IsZero(v);
					var results = new Array(4);
					do {
						var c0 = 0, c1 = 0;
						if (results[0] = CCT.mathRaycastSegment(ls1[0], dir, ls2)) {
							c0 = 1;
							if (!result)
								result = results[0];
						}
						if (results[1] = CCT.mathRaycastSegment(ls1[1], dir, ls2)) {
							c1 = 1;
							if (!result || result.distance > results[1].distance)
								result = results[1];
						}
						if (is_parallel && (c0 || c1))
							break;
						else if (c0 && c1)
							break;
						neg_dir = dir.clone().negate();
						if (results[2] = CCT.mathRaycastSegment(ls2[0], neg_dir, ls1)) {
							if (!result || result.distance > results[2].distance) {
								result = results[2];
								result.hit_point = ls2[0].clone();
							}
						}
						if (results[3] = CCT.mathRaycastSegment(ls2[1], neg_dir, ls1)) {
							if (!result || result.distance > results[3].distance) {
								result = results[3];
								result.hit_point = ls2[1].clone();
							}
						}
					} while (false);
					if (result && is_parallel) {
						var new_ls1 = [
							ls1[0].clone().addScaledVector(dir, result.distance),
							ls1[1].clone().addScaledVector(dir, result.distance)
						];
						res = CCT.overlapSegmentIntersectSegment(new_ls1, ls2);
						if (res && 2 === res.code)
							result.hit_point = undefined;
					}
				}
				return result;
			}
		}
	},
	/**
	 * line segment cast to triangle
	 * @param ls
	 * @param dir
	 * @param tri
	 * @returns {*}
	 */
	mathSegmentcastTriangle : function (ls, dir, tri) {
		var result;
		if (!(result = CCT.mathSegmentcastPlane(ls, dir, tri[0], CCT.mathPlaneNormalByVertices3(tri))))
			return null;
		else if (!result.hit_point) {
			var c = [false, false];
			for (var i = 0; i < 2; ++i) {
				var test_p = ls[i].clone().addScaledVector(dir, result.distance);
				c[i] = (CCT.mathTriangleHasPoint(tri, test_p) !== null);
			}
			if (c[0] && c[1])
				return result;
		}
		else if (CCT.mathTriangleHasPoint(tri, result.hit_point))
			return result;

		var results = new Array(3);
		result = null;
		for (i = 0; i < 3; ++i) {
			var edge = [
				tri[i % 3],
				tri[(i + 1) % 3]
			];
			if (!(results[i] = CCT.mathSegmentcastSegment(ls, dir, edge))) {
				continue;
			}
			if (!result) {
				result = results[i];
			}
			else {
				var cmp = CCT.fcmpf(result.distance, results[i].distance, CCT.EPSILON);
				if (0 === cmp) {
					if (!results[i].hit_point || !result.hit_point ||
						!CCT.mathVec3Equal(result.hit_point, results[i].hit_point))
					{
						result.hit_point = undefined;
					}
					break;
				}
				else if (cmp > 0)
					result = results[i];
			}
		}
		return result;
	},
	/**
	 * line segment cast to sphere
	 * @param ls
	 * @param dir
	 * @param center
	 * @param radius
	 * @returns {*}
	 */
	mathSegmentcastSphere : function (ls, dir, center, radius) {
		var res = CCT.mathSphereIntersectSegment(center, radius, ls);
		if (res && 1 === res.code) {
			return {
				distance : 0.0,
				hit_point : res.p
			};
		}
		else if (res && 2 === res.code) {
			return { distance : 0.0 };
		}
		else {
			var lsdir = new CCT.Vector3().subVectors(ls[1], ls[0]);
			var N = new CCT.Vector3().crossVectors(lsdir, dir);
			var results = new Array(2), result;
			if (CCT.mathVec3IsZero(N)) {
				if (!(results[0] = CCT.mathRaycastSphere(ls[0], dir, center, radius)))
					return null;
				if (!(results[1] = CCT.mathRaycastSphere(ls[1], dir, center, radius)))
					return null;
				return results[0].distance < results[1].distance ? results[0] : results[1];
			}
			else {
				res = CCT.mathSphereIntersectPlane(center, radius, ls[0], N.normalize());
				if (!res)
					return null;
				lsdir.normalize();
				var pj = CCT.mathPointProjectionLine(res.p, ls[0], lsdir);
				var cmp = CCT.fcmpf(pj.v.lengthSq(), res.r * res.r, CCT.EPSILON);
				if (cmp >= 0) {
					var new_ls = new Array(2), d;
					if (cmp > 0) {
						var cos_theta = pj.v.dot(dir);
						if (CCT.fcmpf(cos_theta, 0.0, CCT.EPSILON) <= 0)
							return null;
						d = pj.v.length();
						pj.v.normalize();
						cos_theta = pj.v.dot(dir);
						d -= res.r;
						pj.p.addScaledVector(pj.v, d);
						d /= cos_theta;
						new_ls[0] = ls[0].clone().addScaledVector(dir, d);
						new_ls[1] = ls[1].clone().addScaledVector(dir, d);
					}
					else {
						d = 0.0;
						new_ls[0] = ls[0];
						new_ls[1] = ls[1];
					}
					if (CCT.mathSegmentHasPoint(new_ls, pj.p)) {
						return {
							distance : d,
							hit_point : pj.p,
							hit_normal : pj.v
						};
					}
				}
				result = null;
				for (i = 0; i < 2; ++i) {
					results[i] = CCT.mathRaycastSphere(ls[i], dir, center, radius);
					if (results[i] && (!result || result.distance > results[i].distance))
					{
						result = results[i];
					}
				}
				return result;
			}
		}
	},
	/**
	 *  segment cast to capsule
	 * @param ls
	 * @param dir
	 * @param cp_o
	 * @param cp_axis
	 * @param cp_radius
	 * @param cp_half_height
	 * @returns {*}
	 */
	mathSegmentcastCapsule : function (ls, dir, cp_o, cp_axis, cp_radius, cp_half_height) {
		var i, res = CCT.mathSegmentIntersectCapsule(ls, cp_o, cp_axis, cp_radius, cp_half_height);
		if (res && 1 === res.code) {
			return {
				distance : 0.0,
				hit_point : res.p
			};
		}
		else if (res && 2 === res.code) {
			return {
				distance : 0.0
			};
		}
		else {
			var lsdir = ls[1].clone().sub(ls[0]);
			var N = new CCT.Vector3().crossVectors(lsdir, dir);
			if (CCT.mathVec3IsZero(N)) {
				var result0, result1;
				result0 = CCT.mathRaycastCapsule(ls[0], dir, cp_o, cp_axis, cp_radius, cp_half_height);
				if (!result0)
					return null;
				result1 = CCT.mathRaycastCapsule(ls[1], dir, cp_o, cp_axis, cp_radius, cp_half_height);
				if (!result1)
					return null;
				return result0.distance < result1.distance ? result0 : result1;
			}
			else {
				var result;
				lsdir.normalize();
				res = CCT.mathLineClosestLine(cp_o, cp_axis, ls[0], lsdir);
				if (res.code && CCT.fcmpf(res.min_d, cp_radius, CCT.EPSILON) > 0)
				{
					var closest_p = [
						cp_o.clone().addScaledVector(cp_axis, res.dir_d[0]),
						ls[0].clone().addScaledVector(lsdir, res.dir_d[1])
					];
					var closest_v = closest_p[1].clone().sub(closest_p[0]).normalize();
					var plane_v = closest_p[0].clone().addScaledVector(closest_v, cp_radius);
					result = CCT.mathSegmentcastPlane(ls, dir, plane_v, closest_v);
					if (!result)
						return null;
					else {
						var new_ls = [
							ls[0].clone().addScaledVector(dir, result.distance),
							ls[1].clone().addScaledVector(dir, result.distance)
						];
						res = CCT.mathSegmentIntersectCapsule(new_ls, cp_o, cp_axis, cp_radius, cp_half_height);
						if (res && 1 === res.code) {
							result.hit_point = res.p;
							return result;
						}
						else if (res && 2 === res.code) {
							result.hit_point = null;
							return result;
						}
					}
				}
				result = null;
				for (i = 0; i < 2; ++i) {
					var sphere_o = cp_o.clone().addScaledVector(cp_axis, i ? cp_half_height : -cp_half_height);
					var result_temp = CCT.mathSegmentcastSphere(ls, dir, sphere_o, cp_radius);
					if (!result_temp)
						continue;
					if (!result || result.distance > result_temp.distance)
					{
						result = result_temp;
					}
				}
				for (i = 0; i < 2; ++i) {
					var result_temp = CCT.mathRaycastCapsule(ls[i], dir, cp_o, cp_axis, cp_radius, cp_half_height);
					if (!result_temp)
						continue;
					if (!result || result.distance > result_temp.distance)
					{
						result = result_temp;
					}
				}
				return result;
			}
		}
	},
	/**
	 * triangle cast to plane
	 * @param tri
	 * @param dir
	 * @param vertice
	 * @param normal
	 * @returns {*}
	 */
	mathTrianglecastPlane : function (tri, dir, vertice, normal) {
		var results = new Array(3), result = null;
		for (var i = 0; i < 3; ++i) {
			var ls = [
				tri[i % 3],
				tri[(i + 1) % 3]
			];
			if (results[i] = CCT.mathSegmentcastPlane(ls, dir, vertice, normal)) {
				if (!result)
					result = results[i];
				else {
					var cmp = CCT.fcmpf(result.distance, results[i].distance, CCT.EPSILON);
					if (0 === cmp) {
						if (!results[i].hit_point || !result.hit_point ||
							!CCT.mathVec3Equal(result.hit_point, results[i].hit_point))
						{
							result.hit_point = null;
						}
					}
					else if (cmp > 0)
						result = results[i];
				}
			}
		}
		return result;
	},
	/**
	 * triangle cast to triangle
	 * @param tri1
	 * @param dir
	 * @param tri2
	 * @returns {*}
	 */
	mathTrianglecastTriangle : function (tri1, dir, tri2) {
		var results = new Array(6), result = null;
		var cmp, ls;
		for (var i = 0; i < 3; ++i) {
			ls = [
				tri1[i % 3],
				tri1[(i + 1) % 3]
			];
			if (!(results[i] = CCT.mathSegmentcastTriangle(ls, dir, tri2)))
				continue;
			if (!result) {
				result = results[i];
				continue;
			}
			cmp = CCT.fcmpf(result.distance, results[i].distance, CCT.EPSILON);
			if (cmp > 0) {
				result = results[i];
			}
			else if (0 === cmp && !results[i].hit_point) {
				result = results[i];
			}
		}
		var neg_dir = dir.clone().negate();
		for (; i < 6; ++i) {
			ls = [
				tri2[i % 3],
				tri2[(i + 1) % 3]
			];
			if (!(results[i] = CCT.mathSegmentcastTriangle(ls, neg_dir, tri1)))
				continue;
			if (!result) {
				result = results[i];
				continue;
			}
			cmp = CCT.fcmpf(result.distance, results[i].distance, CCT.EPSILON);
			if (cmp > 0) {
				result = results[i];
			}
			else if (0 === cmp && !results[i].hit_point) {
				result = results[i];
			}
		}
		return result;
	},
	/**
	 * AABB cast to plane
	 * @param o
	 * @param half
	 * @param dir
	 * @param vertice
	 * @param normal
	 * @returns {*}
	 */
	mathAABBcastPlane : function (o, half, dir, vertice, normal) {
		if (CCT.mathAABBIntersectPlane(o, half, vertice, normal)) {
			return { distance : 0.0 };
		}
		else {
			var result = null;
			var v = CCT.AABBVertices(o, half);
			for (var i = 0; i < CCT.Box_Triangle_Vertices_Indices.length; i += 3) {
				var tri = [
					v[CCT.Box_Triangle_Vertices_Indices[i]],
					v[CCT.Box_Triangle_Vertices_Indices[i + 1]],
					v[CCT.Box_Triangle_Vertices_Indices[i + 2]]
				];
				var result_temp = CCT.mathTrianglecastPlane(tri, dir, vertice, normal);
				if (!result_temp)
					continue;
				if (!result || result.distance > result_temp.distance) {
					result = result_temp;
				}
			}
			return result;
		}
	},
	/**
	 * AABB cast to AABB
	 * @param o1
	 * @param half1
	 * @param dir
	 * @param o2
	 * @param half2
	 * @returns {*}
	 */
	mathAABBcastAABB : function (o1, half1, dir, o2, half2) {
		if (CCT.mathAABBIntersectAABB(o1, half1, o2, half2)) {
			return { distance: 0.0 };
		}
		else {
			var result = null;
			var v1 = CCT.AABBVertices(o1, half1);
			var v2 = CCT.AABBVertices(o2, half2);
			for (var i = 0; i < CCT.Box_Triangle_Vertices_Indices.length; i += 3) {
				var tri1 = [
					v1[CCT.Box_Triangle_Vertices_Indices[i]],
					v1[CCT.Box_Triangle_Vertices_Indices[i + 1]],
					v1[CCT.Box_Triangle_Vertices_Indices[i + 2]]
				];
				for (var j = 0; j < CCT.Box_Triangle_Vertices_Indices.length; j += 3) {
					var tri2 = [
						v2[CCT.Box_Triangle_Vertices_Indices[j]],
						v2[CCT.Box_Triangle_Vertices_Indices[j + 1]],
						v2[CCT.Box_Triangle_Vertices_Indices[j + 2]]
					];
					var result_temp = CCT.mathTrianglecastTriangle(tri1, dir, tri2);
					if (!result_temp)
						continue;

					if (!result || result.distance > result_temp.distance) {
						result = result_temp;
					}
				}
			}
			return result;
		}
	},
	/**
	 * sphere cast to plane
	 * @param o
	 * @param radius
	 * @param dir
	 * @param plane_v
	 * @param plane_n
	 * @returns {*}
	 */
	mathSpherecastPlane : function (o, radius, dir, plane_v, plane_n) {
		var pp = CCT.mathPointProjectionPlane(o, plane_v, plane_n);
		var cmp = CCT.fcmpf(pp.d * pp.d, radius * radius, CCT.EPSILON);
		if (cmp < 0) {
			return { distance : 0.0 };
		}
		else if (0 === cmp) {
			return {
				distance : 0.0,
				hit_point : o.clone().addScaledVector(plane_n, pp.d)
			};
		}
		else {
			var cos_theta = plane_n.dot(dir);
			if (CCT.fcmpf(cos_theta, 0.0, CCT.EPSILON) === 0)
				return null;
			var d = pp.d / cos_theta;
			if (CCT.fcmpf(d, 0.0, CCT.EPSILON) < 0)
				return null;
			var pp_d_abs = CCT.fcmpf(pp.d, 0.0, CCT.EPSILON) > 0 ? pp.d : -pp.d;
			d -= radius / pp_d_abs * d;
			var result = { distance : d, hit_normal : plane_n };
			result.hit_point = o.clone().addScaledVector(dir, d);
			if (CCT.fcmpf(pp.d, 0.0, CCT.EPSILON) < 0)
				result.hit_point.addScaledVector(plane_n, -radius);
			else
				result.hit_point.addScaledVector(plane_n, radius);
			return result;
		}
	},
	/**
	 * sphere cast to sphere
	 * @param o1
	 * @param r1
	 * @param dir
	 * @param o2
	 * @param r2
	 * @returns {*}
	 */
	mathSpherecastSphere : function (o1, r1, dir, o2, r2) {
		var result = CCT.mathRaycastSphere(o1, dir, o2, r1 + r2);
		if (!result)
			return null;
		result.hit_normal = result.hit_point.clone().sub(o2).normalize();
		result.hit_point.addScaledVector(result.hit_normal, -r1);
		return result;
	},
	/**
	 * sphere cast to some triangles in same plane
	 * @param o
	 * @param radius
	 * @param dir
	 * @param plane_n
	 * @param vertices
	 * @param indices
	 * @returns {*}
	 */
	mathSpherecastTrianglesPlane : function (o, radius, dir, plane_n, vertices, indices) {
		var res = CCT.mathSphereIntersectTrianglesPlane(o, radius, plane_n, vertices, indices);
		if (res) {
			return { distance: 0.0 };
		}
		var result = CCT.mathSpherecastPlane(o, radius, dir, vertices[indices[0]], plane_n);
		if (result) {
			var i = 0;
			if (result.hit_point) {
				while (i < indices.length) {
					var tri = [
						vertices[indices[i++]],
						vertices[indices[i++]],
						vertices[indices[i++]]
					];
					if (CCT.mathTriangleHasPoint(tri, result.hit_point))
						return result;
				}
			}
			result = null;
			var neg_dir = dir.clone().negate();
			for (i = 0; i < indices.length; i += 3) {
				for (var j = 0; j < 3; ++j) {
					var edge = [
						vertices[indices[j % 3 + i]],
						vertices[indices[(j + 1) % 3 + i]]
					];
					var result_temp = CCT.mathSegmentcastSphere(edge, neg_dir, o, radius);
					if (!result_temp)
						continue;
					if (!result || result.distance > result_temp.distance)
					{
						result = result_temp;
					}
				}
			}
			if (result && result.hit_point) {
				result.hit_point.addScaledVector(dir, result.distance);
			}
		}
		return result;
	},
	/**
	 * sphere cast to AABB
	 * @param o
	 * @param radius
	 * @param dir
	 * @param center
	 * @param half
	 * @returns {*}
	 */
	mathSpherecastAABB : function (o, radius, dir, center, half) {
		if (CCT.mathAABBIntersectSphere(center, half, o, radius)) {
			return { distance : 0.0 };
		}
		else {
			var result = null;
			var v = CCT.AABBVertices(center, half);
			for (var i = 0, j = 0; i < CCT.Box_Triangle_Vertices_Indices.length; i += 6, ++j) {
				var plane_n = new CCT.Vector3().fromArray(CCT.AABB_Plane_Normal[j]);
				var result_temp = CCT.mathSpherecastTrianglesPlane(o, radius, dir,
					plane_n, v, CCT.Box_Triangle_Vertices_Indices.slice(i, i+ 6));
				if (!result_temp)
					continue;
				if (!result || result.distance > result_temp.distance)
				{
					result = result_temp;
				}
			}
			return result;
		}
	},
	/**
	 * sphere cast to capsule
	 * @param sp_o
	 * @param sp_radius
	 * @param dir
	 * @param cp_o
	 * @param cp_axis
	 * @param cp_radius
	 * @param cp_half_height
	 * @returns {*|{distance, hit_point}}
	 */
	mathSpherecastCapsule : function (sp_o, sp_radius, dir, cp_o, cp_axis, cp_radius, cp_half_height) {
		var result = CCT.mathRaycastCapsule(sp_o, dir, cp_o, cp_axis, sp_radius + cp_radius, cp_half_height);
		if (result) {
			var new_sp_o = sp_o.clone().addScaledVector(dir, result.distance);
			var res = CCT.mathSphereIntersectCapsule(new_sp_o, sp_radius, cp_o, cp_axis, cp_radius, cp_half_height);
			if (res)
				result.hit_point = res.p;
		}
		return result;
	},
	/**
	 * capsule cast to plane
	 * @param cp_o
	 * @param cp_axis
	 * @param cp_radius
	 * @param cp_half_height
	 * @param dir
	 * @param plane_v
	 * @param plane_n
	 * @returns {*}
	 */
	mathCapsulecastPlane : function (cp_o, cp_axis, cp_radius, cp_half_height, dir, plane_v, plane_n) {
		var res = CCT.mathCapsuleIntersectPlane(cp_o, cp_axis, cp_radius, cp_half_height, plane_v, plane_n);
		if (res && 2 === res.code) {
			return { distance : 0.0 };
		}
		else if (res && 1 === res.code) {
			return {
				distance : 0.0,
				hit_point : res.p
			};
		}
		else {
			var result = null;
			for (var i = 0; i < 2; ++i) {
				var sphere_o = cp_o.clone().addScaledVector(cp_axis, i ? cp_half_height : -cp_half_height);
				var result_temp = CCT.mathSpherecastPlane(sphere_o, cp_radius, dir, plane_v, plane_n);
				if (!result_temp)
					continue;

				if (!result || result.distance > result_temp.distance)
				{
					result = result_temp;
				}
			}
			if (result) {
				var dot = cp_axis.dot(plane_n);
				if (0 === CCT.fcmpf(dot, 0.0, CCT.EPSILON))
					result.hit_point = undefined;
			}
			return result;
		}
	},
	/**
	 * capsule cast to triangles plane
	 * @param cp_o
	 * @param cp_axis
	 * @param cp_radius
	 * @param cp_half_height
	 * @param dir
	 * @param plane_n
	 * @param vertices
	 * @param indices
	 * @returns {*}
	 */
	mathCapsulecastTrianglesPlane : function (cp_o, cp_axis, cp_radius, cp_half_height, dir, plane_n, vertices, indices) {
		if (CCT.mathCapsuleIntersectTrianglesPlane(cp_o, cp_axis, cp_radius, cp_half_height, plane_n, vertices, indices)) {
			return { distance : 0.0 };
		}
		var result = CCT.mathCapsulecastPlane(cp_o, cp_axis, cp_radius, cp_half_height, dir, vertices[indices[0]], plane_n);
		if (result) {
			var i = 0;
			var p;
			if (result.hit_point) {
				p = result.hit_point;
			}
			else if (CCT.fcmpf(result.distance, 0.0, CCT.EPSILON) > 0) {
				p = cp_o.clone().addScaledVector(dir, result.distance);
				p = CCT.mathPointProjectionPlane(p, vertices[indices[0]], plane_n).p;
			}
			if (p) {
				while (i < indices.length) {
					var tri = [
						vertices[indices[i++]],
						vertices[indices[i++]],
						vertices[indices[i++]]
					];
					if (CCT.mathTriangleHasPoint(tri, p))
						return result;
				}
			}
			var neg_dir = dir.clone().negate();
			result = null;
			for (i = 0; i < indices.length; i += 3) {
				for (var j = 0; j < 3; ++j) {
					var edge = [
						vertices[indices[j % 3 + i]],
						vertices[indices[(j + 1) % 3 + i]]
					];
					var result_temp = CCT.mathSegmentcastCapsule(edge, neg_dir, cp_o, cp_axis, cp_radius, cp_half_height);
					if (!result_temp)
						continue;
					if (!result || result.distance > result_temp.distance)
					{
						result = result_temp;
					}
				}
			}
			if (result && result.hit_point)
				result.hit_point.addScaledVector(dir, result.distance);
		}
		return result;
	},
	/**
	 * capsule cast to AABB
	 * @param cp_o
	 * @param cp_axis
	 * @param cp_radius
	 * @param cp_half_height
	 * @param dir
	 * @param aabb_o
	 * @param aabb_half
	 * @returns {*}
	 */
	mathCapsulecastAABB : function (cp_o, cp_axis, cp_radius, cp_half_height, dir, aabb_o, aabb_half) {
		if (CCT.mathAABBHasPoint(aabb_o, aabb_half, cp_o) ||
			CCT.mathCapsuleHasPoint(cp_o, cp_axis, cp_radius, cp_half_height, aabb_o))
		{
			return { distance : 0.0 };
		}
		else {
			var result = null;
			var v = CCT.AABBVertices(aabb_o, aabb_half);
			for (var i = 0, j = 0; i < CCT.Box_Triangle_Vertices_Indices.length; i += 6, ++j) {
				var plane_n = new CCT.Vector3().fromArray(CCT.AABB_Plane_Normal[j]);
				var result_temp = CCT.mathCapsulecastTrianglesPlane(cp_o, cp_axis, cp_radius, cp_half_height, dir,
					plane_n, v, CCT.Box_Triangle_Vertices_Indices.slice(i, i + 6));
				if (!result_temp)
					continue;
				if (!result || result.distance > result_temp.distance)
				{
					result = result_temp;
				}
			}
			return result;
		}
	},
	/**
	 * capsule cast to capsule
	 * @param cp1_o
	 * @param cp1_axis
	 * @param cp1_radius
	 * @param cp1_half_height
	 * @param dir
	 * @param cp2_o
	 * @param cp2_axis
	 * @param cp2_radius
	 * @param cp2_half_height
	 * @returns {*}
	 */
	mathCapsulecastCapsule : function (cp1_o, cp1_axis, cp1_radius, cp1_half_height, dir, cp2_o, cp2_axis, cp2_radius, cp2_half_height) {
		if (CCT.mathCapsuleIntersectCapsule(cp1_o, cp1_axis, cp1_radius, cp1_half_height, cp2_o, cp2_axis, cp2_radius, cp2_half_height)) {
			return { distance : 0.0 };
		}
		else {
			var sphere_o;
			var N = new CCT.Vector3().crossVectors(cp1_axis, dir);
			if (CCT.mathVec3IsZero(N)) {
				sphere_o = cp1_o.clone().addScaledVector(cp1_axis, cp1_half_height);
				var result0 = CCT.mathSpherecastCapsule(sphere_o, cp1_radius, dir, cp2_o, cp2_axis, cp2_radius, cp2_half_height);
				if (!result0)
					return null;
				sphere_o = cp1_o.clone().addScaledVector(cp1_axis, -cp1_half_height);
				var result1 = CCT.mathSpherecastCapsule(sphere_o, cp1_radius, dir, cp2_o, cp2_axis, cp2_radius, cp2_half_height);
				if (!result1)
					return null;
				return result0.distance < result1.distance ? result0 : result1;
			}
			else {
				var result;
				var i, res = CCT.mathLineClosestLine(cp1_o, cp1_axis, cp2_o, cp2_axis);
				if (res.code && CCT.fcmpf(res.min_d, cp1_radius + cp2_radius, CCT.EPSILON) > 0) {
					var closest_p = [
						cp1_o.clone().addScaledVector(cp1_axis, res.dir_d[0]),
						cp2_o.clone().addScaledVector(cp2_axis, res.dir_d[1])
					];
					var closest_v = closest_p[1].clone().sub(closest_p[0]).normalize();
					var plane_p = cp2_o.clone().addScaledVector(closest_v, -cp2_radius);
					result = CCT.mathCapsulecastPlane(cp1_o, cp1_axis, cp1_radius, cp1_half_height, dir, plane_p, closest_v);
					if (!result)
						return null;
					else {
						var new_cp1_o = cp1_o.clone().addScaledVector(dir, result.distance);
						if (CCT.mathCapsuleIntersectCapsule(new_cp1_o, cp1_axis, cp1_radius, cp1_half_height, cp2_o, cp2_axis, cp2_radius, cp2_half_height)) {
							result.hit_point = undefined;
							return result;
						}
					}
				}
				var neg_dir = dir.clone().negate();
				result = null;
				for (i = 0; i < 2; ++i) {
					sphere_o = cp2_o.clone().addScaledVector(cp2_axis, i ? cp2_half_height : -cp2_half_height);
					var result_temp = CCT.mathSpherecastCapsule(sphere_o, cp2_radius, neg_dir, cp1_o, cp1_axis, cp1_radius, cp1_half_height);
					if (!result_temp)
						continue;
					if (!result || result.distance > result_temp.distance)
						result = result_temp;
				}
				for (i = 0; i < 2; ++i) {
					sphere_o = cp1_o.clone().addScaledVector(cp1_axis, i ? cp1_half_height : -cp1_half_height);
					var result_temp = CCT.mathSpherecastCapsule(sphere_o, cp1_radius, dir, cp2_o, cp2_axis, cp2_radius, cp2_half_height);
					if (!result_temp)
						continue;
					if (!result || result.distance > result_temp.distance) {
						result = result_temp;
					}
				}
				if (result) {
					result.hit_point = undefined;
				}
				return result;
			}
		}
	}
};

CCT.Vector3.prototype.clone = function () {
	return new CCT.Vector3(this.x, this.y, this.z);
};
CCT.Vector3.prototype.addScaledVector = function (v, s) {
	this.x += v.x * s;
	this.y += v.y * s;
	this.z += v.z * s;
	return this;
};
CCT.Vector3.prototype.sub = function (v) {
	this.x -= v.x;
	this.y -= v.y;
	this.z -= v.z;
	return this;
};
CCT.Vector3.prototype.subVectors = function (a, b) {
	this.x = a.x - b.x;
	this.y = a.y - b.y;
	this.z = a.z - b.z;
	return this;
};
CCT.Vector3.prototype.negate = function () {
	this.x = -this.x;
	this.y = -this.y;
	this.z = -this.z;
	return this;
};
CCT.Vector3.prototype.dot = function (v) {
	return this.x * v.x + this.y * v.y + this.z * v.z;
};
CCT.Vector3.prototype.lengthSq = function () {
	return this.x * this.x + this.y * this.y + this.z * this.z;
};
CCT.Vector3.prototype.length = function () {
	return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
};
CCT.Vector3.prototype.multiplyScalar = function (scalar) {
	this.x *= scalar;
	this.y *= scalar;
	this.z *= scalar;
	return this;
};
CCT.Vector3.prototype.divideScalar = function (scalar) {
	return this.multiplyScalar(1 / scalar);
};
CCT.Vector3.prototype.normalize = function () {
	return this.divideScalar(this.length() || 1);
};
CCT.Vector3.prototype.crossVectors = function (a, b) {
	var ax = a.x, ay = a.y, az = a.z;
	var bx = b.x, by = b.y, bz = b.z;
	this.x = ay * bz - az * by;
	this.y = az * bx - ax * bz;
	this.z = ax * by - ay * bx;
	return this;
};
CCT.Vector3.prototype.fromArray = function (array, offset) {
	if (offset === undefined)
		offset = 0;
	this.x = array[offset];
	this.y = array[offset + 1];
	this.z = array[offset + 2];
	return this;
};
