#include <cmath>

namespace sMath {
	struct float3d;
	struct matrix4x4;
	matrix4x4 ProjectionMatrix(float fNear, float fFar, float fFov, int screenWidth, int screenHeight);
	matrix4x4 PointAtMatrix(float3d& pos, float3d& target, float3d& up);
	matrix4x4 PointAtMatrixInv(matrix4x4);
	matrix4x4 xRotateMat(float);
	matrix4x4 yRotateMat(float);
	matrix4x4 zRotateMat(float);
	float3d vecMultMat4x4(float3d, matrix4x4);
	float vecLength(float3d);
	float dotProduct(float3d, float3d);
	float vecToVecAngleR(float3d, float3d);
	float vecToVecAngleD(float3d, float3d);
	float3d VecIntersectPlane(float3d& plane_p, float3d& plane_n, float3d& lineStart, float3d& lineEnd);
	float3d buildVec(float3d, float3d);
	float3d buildNormal(float3d, float3d, float3d, float3d);
	float3d vecNormalize(float3d);
	float3d crossProduct(float3d, float3d);
	float3d scalMultVec(float, float3d);
	float3d vecPlusVec(float3d, float3d);
	float3d xUnit();
	float3d yUnit();
	float3d zUnit();
}




	// Point/Vector structure
	struct sMath::float3d
	{
		float x;
		float y;
		float z;
		float w = 1.0f;
		float3d() {};
		float3d(float a, float b, float c) :x(a), y(b), z(c) {};
	};

	// Matrix structure
	struct sMath::matrix4x4 
	{
		float m[4][4] = { 0 };
	};

	// Returns the projection matrix
	sMath::matrix4x4 sMath::ProjectionMatrix(float fNear, float fFar, float fFov, int nScreenWidth, int nScreenHeight)
	{
		sMath::matrix4x4 mat;
		float fFovRad = 1.0f / tanf(fFov * 0.5f / 180.0f * 3.14159f);
		float fAspectRatio = (float)nScreenHeight / (float)nScreenWidth;
		mat.m[0][0] = fAspectRatio * fFovRad;
		mat.m[1][1] = fFovRad;
		mat.m[2][2] = fFar / (fFar - fNear);
		mat.m[3][2] = -fFar * fNear / (fFar - fNear);
		mat.m[2][3] = 1.0f;
		return mat;
	}

	// Returns the point-at matrix
	sMath::matrix4x4 sMath::PointAtMatrix(sMath::float3d& pos, sMath::float3d& target, sMath::float3d& up)
	{
		sMath::matrix4x4 mat;
		sMath::float3d newForward = sMath::buildVec(pos, target);
		sMath::vecNormalize(newForward);

		sMath::float3d a = sMath::scalMultVec(sMath::dotProduct(up, newForward), newForward);
		sMath::float3d newUp = sMath::vecPlusVec(a, sMath::scalMultVec(-1.0f, up));
		newUp = sMath::vecNormalize(newUp);

		sMath::float3d newRight = sMath::crossProduct(newForward, newUp);

		mat.m[0][0] = newRight.x;	mat.m[0][1] = newRight.y;	mat.m[0][2] = newRight.z;	mat.m[0][3] = 0.0f;
		mat.m[1][0] = newUp.x;		mat.m[1][1] = newUp.y;		mat.m[1][2] = newUp.z;		mat.m[1][3] = 0.0f;
		mat.m[2][0] = newForward.x;	mat.m[2][1] = newForward.y;	mat.m[2][2] = newForward.z;	mat.m[2][3] = 0.0f;
		mat.m[3][0] = pos.x;		mat.m[3][1] = pos.y;		mat.m[3][2] = pos.z;		mat.m[3][3] = 1.0f;
		return mat;
	}

	// Returns the inverse of the point-at matrix
	sMath::matrix4x4 sMath::PointAtMatrixInv(sMath::matrix4x4 m)
	{
		float dpTA = m.m[0][0] * m.m[3][0] + m.m[0][1] * m.m[3][1] + m.m[0][2] * m.m[3][2];
		float dpTB = m.m[1][0] * m.m[3][0] + m.m[1][1] * m.m[3][1] + m.m[1][2] * m.m[3][2];
		float dpTC = m.m[2][0] * m.m[3][0] + m.m[2][1] * m.m[3][1] + m.m[2][2] * m.m[3][2];


		sMath::matrix4x4 mat;
		mat.m[0][0] = m.m[0][0];	mat.m[0][1] = m.m[1][0];	mat.m[0][2] = m.m[2][0];	mat.m[0][3] = 0.0f;
		mat.m[1][0] = m.m[0][1];	mat.m[1][1] = m.m[1][1];	mat.m[1][2] = m.m[2][1];	mat.m[1][3] = 0.0f;
		mat.m[2][0] = m.m[0][2];	mat.m[2][1] = m.m[1][2];	mat.m[2][2] = m.m[2][2];	mat.m[2][3] = 0.0f;
		mat.m[3][0] = -dpTA;		mat.m[3][1] = -dpTB;		mat.m[3][2] = -dpTC;		mat.m[3][3] = 1.0f;
		return mat;
	}

	// Multiplies two 4x4 matrices
	sMath::float3d sMath::vecMultMat4x4(sMath::float3d vec, sMath::matrix4x4 mat)
	{
		float3d res;
		res.x = vec.x * mat.m[0][0] + vec.y * mat.m[1][0] + vec.z * mat.m[2][0] + vec.w * mat.m[3][0];
		res.y = vec.x * mat.m[0][1] + vec.y * mat.m[1][1] + vec.z * mat.m[2][1] + vec.w * mat.m[3][1];
		res.z = vec.x * mat.m[0][2] + vec.y * mat.m[1][2] + vec.z * mat.m[2][2] + vec.w * mat.m[3][2];
		res.w = vec.x * mat.m[0][3] + vec.y * mat.m[1][3] + vec.z * mat.m[2][3] + vec.w * mat.m[3][3];

		return res;
	}

	// Returns the rotation matrix (about x axis)
	sMath::matrix4x4 sMath::xRotateMat(float theta)
	{
		sMath::matrix4x4 mat;
		mat.m[0][0] = 1.0f;
		mat.m[1][1] = cosf(theta);
		mat.m[1][2] = sinf(theta);
		mat.m[2][1] = -sinf(theta);
		mat.m[2][2] = cos(theta);
		mat.m[3][3] = 1.0f;
		return mat;
	}

	// Returns the rotation matrix (about y axis)
	sMath::matrix4x4 sMath::yRotateMat(float theta)
	{
		sMath::matrix4x4 mat;
		mat.m[0][0] = cosf(theta);
		mat.m[1][1] = 1.0f;
		mat.m[0][2] = -sinf(theta);
		mat.m[2][0] = sinf(theta);
		mat.m[2][2] = cos(theta);
		mat.m[3][3] = 1.0f;
		return mat;
	}

	// Returns the rotation matrix (about z axis)
	sMath::matrix4x4 sMath::zRotateMat(float theta)
	{
		sMath::matrix4x4 mat;
		mat.m[0][0] = cosf(theta*0.5f);
		mat.m[1][1] = cosf(theta*0.5f);
		mat.m[0][1] = sinf(theta*0.5f);
		mat.m[1][0] = -sinf(theta*0.5f);
		mat.m[2][2] = 1.0f;
		mat.m[3][3] = 1.0f;
		return mat;
	}

	// Computes the lenth of a vector
	float sMath::vecLength(sMath::float3d v)
	{
		return (float)(std::sqrtf(std::powf(v.x, 2.0f) + std::powf(v.y, 2.0f) + std::powf(v.z, 2.0f)));
	}

	// Computes the dot product of two vectors
	float sMath::dotProduct(sMath::float3d a, sMath::float3d b)
	{
		return (float)(a.x * b.x + a.y * b.y + a.z * b.z);
	}

	// Computes the cross product of two vectors
	sMath::float3d sMath::crossProduct(sMath::float3d a, sMath::float3d b)
	{
		sMath::float3d res;
		res.x = (float)a.y * b.z - (float)a.z * b.y;
		res.y = (float)a.z * b.x - (float)a.x * b.z;
		res.z = (float)a.x * b.y - (float)a.y * b.x;
		return res;
	}

	// Returns the angle between two vectors (in radians)
	float sMath::vecToVecAngleR(sMath::float3d a, sMath::float3d b)
	{
		return std::asinf(sMath::vecLength(sMath::crossProduct(a, b)) / (sMath::vecLength(a) * sMath::vecLength(b)));
	}

	// Returns the angle between two vectors (in degrees)
	float sMath::vecToVecAngleD(sMath::float3d a, sMath::float3d b)
	{
		return (std::asinf(vecLength(crossProduct(a, b)) / (vecLength(a) * vecLength(b))) * 180.0f / 3.14159265f);
	}

	// Takes two points and returns a vector built out of them
	sMath::float3d sMath::buildVec(float3d const p1, float3d const p2)
	{
		sMath::float3d res;
		res.x = p2.x - p1.x;
		res.y = p2.y - p1.y;
		res.z = p2.z - p1.z;
		return res;
	}

	//Vector and Plane intersection
	sMath::float3d sMath::VecIntersectPlane(sMath::float3d& plane_p, sMath::float3d& plane_n, sMath::float3d& lineStart, sMath::float3d& lineEnd)
	{

		plane_n = sMath::vecNormalize(plane_n);
		float plane_d = -sMath::dotProduct(plane_n, plane_p);
		float ad = sMath::dotProduct(lineStart, plane_n);
		float bd = sMath::dotProduct(lineEnd, plane_n);
		float t = (-plane_d - ad) / (bd - ad);
		sMath::float3d lineStartToEnd = sMath::vecPlusVec(sMath::scalMultVec(-1.0f, lineStart), lineEnd );
		sMath::float3d lineToIntersect = sMath::scalMultVec(t, lineStartToEnd);
		return sMath::vecPlusVec(lineStart, lineToIntersect);
	}

	//B
uilds a normal to a plane out of four points
	sMath::float3d sMath::buildNormal(sMath::float3d const p1, sMath::float3d const p2, sMath::float3d const p3, sMath::float3d const p4)
	{
		sMath::float3d vec1, vec2, res;
		vec1 = sMath::buildVec(p1, p2);
		vec2 = sMath::buildVec(p3, p4);
		res = sMath::crossProduct(vec1, vec2);
		res = sMath::vecNormalize(res);
		return res;
	}

	// Normalizes a vector
	sMath::float3d sMath::vecNormalize(sMath::float3d v)
	{
		sMath::float3d n;
		float vlen = sMath::vecLength(v);
		n.x = v.x / vlen;
		n.y = v.y / vlen;
		n.z = v.z / vlen;
		return n;
	}

	//Computes vector multiplied by a number
	sMath::float3d sMath::scalMultVec(float n, sMath::float3d v)
	{
		v.x *= n;
		v.y *= n;
		v.z *= n;
		return v;
	}

	//Computes the sum of two vectors
	sMath::float3d sMath::vecPlusVec(sMath::float3d a, sMath::float3d b)
	{
		sMath::float3d res;
		res.x = a.x + b.x;
		res.y = a.y + b.y;
		res.z = a.z + b.z;
		return res;
	}

	//Returns the unit vector of x axis
	sMath::float3d sMath::xUnit()
	{
		sMath::float3d x(1.0f, 0.0f, 0.0f);
		return x;
	}

	//Returns the unit vector of y axis
	sMath::float3d sMath::yUnit()
	{
		sMath::float3d y(0.0f, 1.0f, 0.0f);
		return y;
	}

	//Returns the unit vector of z axis
	sMath::float3d sMath::zUnit()
	{
		sMath::float3d z(0.0f, 0.0f, 1.0f);
		return z;
	}
