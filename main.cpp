#include "olcConsoleGameEngine.h";
#include "sMath.h";
#include <vector>;
#include <fstream>;
#include <strstream>;
#include <string>;
#include <algorithm>;

struct Triangle
{
	sMath::float3d p[3];
	wchar_t sym;
	short col;
};

struct Mesh
{
	std::vector<Triangle> tris;

	bool loadFromObjFile(std::string sFilename)
	{
		std::ifstream f(sFilename);

		if (!f.is_open())
		{
			std::cout << "Error";
			return false;
		}

		std::vector<sMath::float3d> verts;
		while (!f.eof())
		{
			char line[256];
			f.getline(line, 256);
			
			std::strstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				sMath::float3d v;
				s >> junk >> v.x >> v.y >> v.z;
				verts.push_back(v);
			}

			if (line[0] == 'f')
			{
				int f[3];
				/*std::string tmp1, tmp2, tmp3, tmp;
				int i = 0;
				s >> junk >> tmp1 >> tmp2 >> tmp3;

				while (tmp1[i] != '/')
				{
					tmp.push_back(tmp1[i]);
					i++;
				}
				f[0] = stoi(tmp);
				i = 0;
				tmp.clear();

				while (tmp2[i] != '/')
				{
					tmp.push_back(tmp2[i]);
					i++;
				}
				f[1] = stoi(tmp);
				i = 0;
				tmp.clear();

				while (tmp3[i] != '/')
				{
					tmp.push_back(tmp3[i]);
					i++;
				}
				f[2] = stoi(tmp);
				i = 0;
				tmp.clear();*/
				
				s >> junk >> f[0] >> f[1] >> f[2];

				Triangle t;
				t.p[0] = verts[f[0] - 1];

				t.p[1] = verts[f[1] - 1];

				t.p[2] = verts[f[2] - 1];

				tris.push_back(t);
			}
		}

		//float tmpx = 0.0f - tris[0].p[0].x;
		//float tmpy = 0.0f - tris[0].p[0].y;
		//float tmpz = 0.0f - tris[0].p[0].z;

		//for (auto v : tris)
		//{
		//	v.p[0].x += tmpx; v.p[1].x += tmpx; v.p[2].x += tmpx;
		//	v.p[0].y += tmpy; v.p[1].y += tmpy;	v.p[2].y += tmpy;
		//	v.p[0].z += tmpz; v.p[1].z += tmpz;	v.p[2].z += tmpz;
		//}

		return true;
	}
};

int Triangle_ClipAgainstPlane(sMath::float3d plane_p, sMath::float3d plane_n, Triangle& in_tri, Triangle& out_tri1, Triangle& out_tri2)
{
	// Make sure plane normal is indeed normal
	plane_n = sMath::vecNormalize(plane_n);

	// Return signed shortest distance from point to plane, plane normal must be normalised
	auto dist = [&](sMath::float3d& p)
	{
		sMath::float3d n = sMath::vecNormalize(p);
		return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - sMath::dotProduct(plane_n, plane_p));
	};

	// Create two temporary storage arrays to classify points on either side of a plane
	// If distance sign is positive, point lies on the "inside" of a plane
	sMath::float3d* inside_points[3];  int nInsidePointCount = 0;
	sMath::float3d* outside_points[3]; int nOutsidePointCount = 0;

	// Get signed distance of each point in triangle to plane
	float d0 = dist(in_tri.p[0]);
	float d1 = dist(in_tri.p[1]);
	float d2 = dist(in_tri.p[2]);

	if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; }
	if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[1]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[1]; }
	if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[2]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[2]; }

	// Now classify triangle points, and break the input triangle into 
	// smaller output triangles if required. There are four possible
	// outcomes...

	if (nInsidePointCount == 0)
	{
		// All points lie on the outside of plane, so clip whole triangle
		// It ceases to exist

		return 0; // No returned triangles are valid
	}

	if (nInsidePointCount == 3)
	{
		// All points lie on the inside of plane, so do nothing
		// and allow the triangle to simply pass through
		out_tri1 = in_tri;

		return 1; // Just the one returned original triangle is valid
	}

	if (nInsidePointCount == 1 && nOutsidePointCount == 2)
	{
		// Triangle should be clipped. As two points lie outside
		// the plane, the triangle simply becomes a smaller triangle

		// Copy appearance info to new triangle
		out_tri1.col = FG_BLUE;//in_tri.col;
		out_tri1.sym = in_tri.sym;

		// The inside point is valid, so keep that...
		out_tri1.p[0] = *inside_points[0];

		// but the two new points are at the locations where the 
		// original sides of the triangle (lines) intersect with the plane
		out_tri1.p[1] = sMath::VecIntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
		out_tri1.p[2] = sMath::VecIntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1]);

		return 1; // Return the newly formed single triangle
	}
	if (nInsidePointCount == 2 && nOutsidePointCount == 1)
	{
		// Triangle should be clipped. As two points lie inside the plane,
		// the clipped triangle becomes a "quad". Fortunately, we can
		// represent a quad with two new triangles

		// Copy appearance info to new triangles
		out_tri1.col = FG_GREEN;//in_tri.col;
		out_tri1.sym = in_tri.sym;

		out_tri2.col = FG_RED;//in_tri.col;
		out_tri2.sym = in_tri.sym;

		// The first triangle consists of the two inside points and a new
		// point determined by the location where one side of the triangle
		// intersects with the plane
		out_tri1.p[0] = *inside_points[0];
		out_tri1.p[1] = *inside_points[1];
		out_tri1.p[2] = sMath::VecIntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);

		// The second triangle is composed of one of he inside points, a
		// new point determined by the intersection of the other side of the 
		// triangle and the plane, and the newly created point above
		out_tri2.p[0] = *inside_points[1];
		out_tri2.p[1] = out_tri1.p[2];
		out_tri2.p[2] = sMath::VecIntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0]);

		return 2; // Return two newly formed triangles which form a quad
	}
}

CHAR_INFO GetColour(float luminance)
{
	wchar_t sym;
	short bg_col, fg_col;
	int pixel_bw = (int)(luminance * 13.0f);
	switch (pixel_bw)
	{
	case 0: bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID; break;

	case 1: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_QUARTER; break;
	case 2: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_HALF; break;
	case 3: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_THREEQUARTERS; break;
	case 4: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_SOLID; break;

	case 5: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_QUARTER; break;
	case 6: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_HALF; break;
	case 7: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_THREEQUARTERS; break;
	case 8: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_SOLID; break;

	case 9:  bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_QUARTER; break;
	case 10: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_HALF; break;
	case 11: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_THREEQUARTERS; break;
	case 12: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_SOLID; break;
	default:
		bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID;
	}
	CHAR_INFO c;
	c.Attributes = bg_col | fg_col;
	c.Char.UnicodeChar = sym;
	return c;
}
 
class RenderEngine : public olcConsoleGameEngine
{
public:
	RenderEngine()
	{
		m_sAppName = L"3d Engine";
	}

	bool olcConsoleGameEngine::OnUserCreate()
	{
		//SOUTH
		Triangle t1,t2;
		t1.p[0].x = 0.0f; t1.p[1].x = 0.0f; t1.p[2].x = 1.0f;
		t1.p[0].y = 0.0f; t1.p[1].y = 1.0f; t1.p[2].y = 1.0f;
		t1.p[0].z = 0.0f; t1.p[1].z = 0.0f; t1.p[2].z = 0.0f;

		t2.p[0].x = 0.0f; t2.p[1].x = 1.0f; t2.p[2].x = 1.0f;
		t2.p[0].y = 0.0f; t2.p[1].y = 1.0f; t2.p[2].y = 0.0f;
		t2.p[0].z = 0.0f; t2.p[1].z = 0.0f; t2.p[2].z = 0.0f;

		//EAST
		Triangle t3, t4;
		t3.p[0].x = 1.0f; t3.p[1].x = 1.0f; t3.p[2].x = 1.0f;
		t3.p[0].y = 0.0f; t3.p[1].y = 1.0f; t3.p[2].y = 1.0f;
		t3.p[0].z = 0.0f; t3.p[1].z = 0.0f; t3.p[2].z = 1.0f;

		t4.p[0].x = 1.0f; t4.p[1].x = 1.0f; t4.p[2].x = 1.0f;
		t4.p[0].y = 0.0f; t4.p[1].y = 1.0f; t4.p[2].y = 0.0f;
		t4.p[0].z = 0.0f; t4.p[1].z = 1.0f; t4.p[2].z = 1.0f;

		//NORTH
		Triangle t5, t6;
		t5.p[0].x = 1.0f; t5.p[1].x = 1.0f; t5.p[2].x = 0.0f;
		t5.p[0].y = 0.0f; t5.p[1].y = 1.0f; t5.p[2].y = 1.0f;
		t5.p[0].z = 1.0f; t5.p[1].z = 1.0f; t5.p[2].z = 1.0f;

		t6.p[0].x = 1.0f; t6.p[1].x = 0.0f; t6.p[2].x = 0.0f;
		t6.p[0].y = 0.0f; t6.p[1].y = 1.0f; t6.p[2].y = 0.0f;
		t6.p[0].z = 1.0f; t6.p[1].z = 1.0f; t6.p[2].z = 1.0f;

		//WEST
		Triangle t7, t8;
		t7.p[0].x = 0.0f; t7.p[1].x = 0.0f; t7.p[2].x = 0.0f;
		t7.p[0].y = 0.0f; t7.p[1].y = 1.0f; t7.p[2].y = 1.0f;
		t7.p[0].z = 1.0f; t7.p[1].z = 1.0f; t7.p[2].z = 0.0f;

		t8.p[0].x = 0.0f; t8.p[1].x = 0.0f; t8.p[2].x = 0.0f;
		t8.p[0].y = 0.0f; t8.p[1].y = 1.0f; t8.p[2].y = 0.0f;
		t8.p[0].z = 1.0f; t8.p[1].z = 0.0f; t8.p[2].z = 0.0f;

		//TOP
		Triangle t9, t10;
		t9.p[0].x = 0.0f; t9.p[1].x = 0.0f; t9.p[2].x = 1.0f;
		t9.p[0].y = 1.0f; t9.p[1].y = 1.0f; t9.p[2].y = 1.0f;
		t9.p[0].z = 0.0f; t9.p[1].z = 1.0f; t9.p[2].z = 1.0f;

		t10.p[0].x = 0.0f; t10.p[1].x = 1.0f; t10.p[2].x = 1.0f;
		t10.p[0].y = 1.0f; t10.p[1].y = 1.0f; t10.p[2].y = 1.0f;
		t10.p[0].z = 0.0f; t10.p[1].z = 1.0f; t10.p[2].z = 0.0f;

		//BOTTOM
		Triangle t11, t12;
		t11.p[0].x = 1.0f; t11.p[1].x = 0.0f; t11.p[2].x = 0.0f;
		t11.p[0].y = 0.0f; t11.p[1].y = 0.0f; t11.p[2].y = 0.0f;
		t11.p[0].z = 1.0f; t11.p[1].z = 1.0f; t11.p[2].z = 0.0f;

		t12.p[0].x = 1.0f; t12.p[1].x = 0.0f; t12.p[2].x = 1.0f;
		t12.p[0].y = 0.0f; t12.p[1].y = 0.0f; t12.p[2].y = 0.0f;
		t12.p[0].z = 1.0f; t12.p[1].z = 0.0f; t12.p[2].z = 0.0f;

		bool err = CubeMesh.loadFromObjFile("ship.obj");
		if(!err)
			CubeMesh.tris = { t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12 };

		float theta = 0.0f;
		
		ProjMatrix = sMath::ProjectionMatrix(0.1f, 1000.0f, 90.0f, ScreenWidth(), ScreenHeight());
		light_direction.z = -1.0f;
		return true;
	}

	bool olcConsoleGameEngine::OnUserUpdate(float fElapsedTime)
	{
		if (GetKey(VK_UP).bHeld)
		{
			Camera.y += 8.0f * fElapsedTime;
		}

		if (GetKey(VK_DOWN).bHeld)
		{
			Camera.y -= 8.0f * fElapsedTime;
		}

		if (GetKey(VK_LEFT).bHeld)
		{
			Camera.x -= 8.0f * fElapsedTime;
		}

		if (GetKey(VK_RIGHT).bHeld)
		{
			Camera.x += 8.0f * fElapsedTime;
		}

		sMath::float3d Forward = sMath::scalMultVec(8.0f * fElapsedTime, LookDir);

		if (GetKey(L'W').bHeld)
		{
			Camera = sMath::vecPlusVec(Camera, Forward);
		}

		if (GetKey(L'S').bHeld)
		{
			Camera = sMath::vecPlusVec(Camera, sMath::scalMultVec(-1.0f, Forward));
		}

		if (GetKey(L'A').bHeld)
		{
			fYaw -= 2.0f * fElapsedTime;
		}

		if (GetKey(L'D').bHeld)
		{
			fYaw += 2.0f * fElapsedTime;
		}


		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);


		theta += 0.4f * fElapsedTime;
		if (theta > 10000.0f) theta = 0.1f;

		sMath::float3d Up;
		Up.x = 0.0f; Up.y = 1.0f; Up.z = 0.0f;
		sMath::float3d Target;
		Target.x = 0.0f; Target.y = 0.0f; Target.z = 1.0f;
		sMath::matrix4x4 matCameraRot = sMath::yRotateMat(fYaw);

		LookDir = sMath::vecMultMat4x4(Target, matCameraRot);
		Target = sMath::vecPlusVec(Camera, LookDir);


		sMath::matrix4x4 CameraMat = sMath::PointAtMatrix(Camera, Target, Up);
		sMath::matrix4x4 ViewMat = sMath::PointAtMatrixInv(CameraMat);

		std::vector<Triangle> vecTrianglesToRaster;

		for (auto tri : CubeMesh.tris)
		{
			Triangle ProjectedTri,TranslatedTri,zRotatedTri,xRotatedTri,xzRotatedTri,ViewedTriangle;

			zRotatedTri.p[0] = sMath::vecMultMat4x4(tri.p[0], sMath::zRotateMat(theta));
			zRotatedTri.p[1] = sMath::vecMultMat4x4(tri.p[1], sMath::zRotateMat(theta));
			zRotatedTri.p[2] = sMath::vecMultMat4x4(tri.p[2], sMath::zRotateMat(theta));

			///*zRotatedTri.p[0] = sMath::vecMultMat4x4(zRotatedTri.p[0], sMath::yRotateMat(0.6 * theta));
			//zRotatedTri.p[1] = sMath::vecMultMat4x4(zRotatedTri.p[1], sMath::yRotateMat(0.6 * theta));
			//zRotatedTri.p[2] = sMath::vecMultMat4x4(zRotatedTri.p[2], sMath::yRotateMat(0.6 * theta));*/

			xzRotatedTri.p[0] = sMath::vecMultMat4x4(zRotatedTri.p[0], sMath::xRotateMat(0.5f * theta));
			xzRotatedTri.p[1] = sMath::vecMultMat4x4(zRotatedTri.p[1], sMath::xRotateMat(0.5f * theta));
			xzRotatedTri.p[2] = sMath::vecMultMat4x4(zRotatedTri.p[2], sMath::xRotateMat(0.5f * theta));

			TranslatedTri = xzRotatedTri;
			TranslatedTri.p[0].z += 10.0f;
			TranslatedTri.p[1].z += 10.0f;
			TranslatedTri.p[2].z += 10.0f;


			sMath::float3d normal = sMath::buildNormal(TranslatedTri.p[0], TranslatedTri.p[1], TranslatedTri.p[0], TranslatedTri.p[2]);
			if (sMath::dotProduct(normal,sMath::buildVec(Camera,TranslatedTri.p[0])) < 0.0f )
			{
				light_direction = sMath::vecNormalize(light_direction);
				float dp = sMath::dotProduct(normal, light_direction);
				CHAR_INFO c = GetColour(dp);

				TranslatedTri.col = c.Attributes;
				TranslatedTri.sym = c.Char.UnicodeChar;

				ViewedTriangle.p[0] = sMath::vecMultMat4x4(TranslatedTri.p[0], ViewMat);
				ViewedTriangle.p[1] = sMath::vecMultMat4x4(TranslatedTri.p[1], ViewMat);
				ViewedTriangle.p[2] = sMath::vecMultMat4x4(TranslatedTri.p[2], ViewMat);
				ViewedTriangle.sym = TranslatedTri.sym;
				ViewedTriangle.col = TranslatedTri.col;

				int nClippedTriangles = 0;
				Triangle clipped[2];
				nClippedTriangles = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, ViewedTriangle, clipped[0], clipped[1]);

			

				for(int n = 0; n < nClippedTriangles; n++)
				{
					ProjectedTri.p[0] = sMath::vecMultMat4x4(clipped[n].p[0], ProjMatrix);
					ProjectedTri.p[1] = sMath::vecMultMat4x4(clipped[n].p[1], ProjMatrix);
					ProjectedTri.p[2] = sMath::vecMultMat4x4(clipped[n].p[2], ProjMatrix);
					ProjectedTri.col = clipped[n].col;
					ProjectedTri.sym = clipped[n].sym;

					ProjectedTri.p[0].x /= ProjectedTri.p[0].w; ProjectedTri.p[1].x /= ProjectedTri.p[1].w; ProjectedTri.p[2].x /= ProjectedTri.p[2].w;
					ProjectedTri.p[0].y /= ProjectedTri.p[0].w;	ProjectedTri.p[1].y /= ProjectedTri.p[1].w;	ProjectedTri.p[2].y /= ProjectedTri.p[2].w;
					ProjectedTri.p[0].z /= ProjectedTri.p[0].w; ProjectedTri.p[1].z /= ProjectedTri.p[1].w;	ProjectedTri.p[2].z /= ProjectedTri.p[2].w;

					ProjectedTri.p[0].x += 1.0f; ProjectedTri.p[0].y += 1.0f;
					ProjectedTri.p[1].x += 1.0f; ProjectedTri.p[1].y += 1.0f;
					ProjectedTri.p[2].x += 1.0f; ProjectedTri.p[2].y += 1.0f;

					ProjectedTri.p[0].x *= 0.5f * (float)ScreenWidth();
					ProjectedTri.p[0].y *= 0.5f * (float)ScreenHeight();
					ProjectedTri.p[1].x *= 0.5f * (float)ScreenWidth();
					ProjectedTri.p[1].y *= 0.5f * (float)ScreenHeight();
					ProjectedTri.p[2].x *= 0.5f * (float)ScreenWidth();
					ProjectedTri.p[2].y *= 0.5f * (float)ScreenHeight();

					ProjectedTri.col = TranslatedTri.col;
					ProjectedTri.sym = TranslatedTri.sym;
				



					vecTrianglesToRaster.push_back(ProjectedTri);				
				}
			}
		}

		std::sort(vecTrianglesToRaster.begin(), vecTrianglesToRaster.end(), [](Triangle& t1, Triangle& t2)
			{
				float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
				float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
				return z1 > z2;
			});

		for (auto& triToRaster : vecTrianglesToRaster)
		{

			Triangle clipped[2];
			std::list<Triangle> listTriangles;
			listTriangles.push_back(triToRaster);
			int nNewTriangles = 1;

			for (int p = 0; p < 4; p++)
			{
				int nTrisToAdd = 0;
				while (nNewTriangles > 0)
				{
					// Take triangle from front of queue
					Triangle test = listTriangles.front();
					listTriangles.pop_front();
					nNewTriangles--;


					switch (p)
					{
					case 0:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 1:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, (float)ScreenHeight() - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 2:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 3:	nTrisToAdd = Triangle_ClipAgainstPlane({ (float)ScreenWidth() - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					}

					// Clipping may yield a variable number of triangles, so
					// add these new ones to the back of the queue for subsequent
					// clipping against next planes
					for (int w = 0; w < nTrisToAdd; w++)
						listTriangles.push_back(clipped[w]);
				}
				nNewTriangles = listTriangles.size();
			}

			for (auto& t : listTriangles)
			{
				FillTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, t.sym, t.col);
				DrawTriangle(t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, PIXEL_SOLID, FG_BLACK);
			}


		}


		return true;
	}
private:
	Mesh CubeMesh;
	sMath::matrix4x4 ProjMatrix;
	float theta = 0.0f;
	sMath::float3d Camera;
	sMath::float3d LookDir;
	sMath::float3d light_direction;
	float fYaw;
};

int main()
{
	RenderEngine engine;

	if (engine.ConstructConsole(256, 240, 2, 2))
		engine.Start();
	return 0;
}
