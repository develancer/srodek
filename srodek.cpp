 ///////////////////////////////////////////////////////////////////////
// various include declarations ///////////////////////////////////////

#include <cstdio>
#include <exception>
#include <mpi.h>

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

#include "shapefil.h"

#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Gnomonic.hpp>
#include <GeographicLib/PolygonArea.hpp>
#include <GeographicLib/TransverseMercator.hpp>
#include <GeographicLib/AlbersEqualArea.hpp>

using GeographicLib::DMS;
using GeographicLib::Geocentric;
using GeographicLib::Geodesic;
using GeographicLib::PolygonArea;

 ///////////////////////////////////////////////////////////////////////
// internal type definitions //////////////////////////////////////////

typedef boost::geometry::model::d2::point_xy<double> Point;
typedef boost::geometry::model::polygon<Point> Polygon;
typedef boost::geometry::model::box<Point> Rectangle;
typedef std::unique_ptr<SHPObject, decltype(&SHPDestroyObject)> ShapeObjectHandle;

 ///////////////////////////////////////////////////////////////////////
// geographic constants for data set //////////////////////////////////

const double GRS80_A = 6378137.0;
const double GRS80_F = 1.0/298.257222101;

const double MERCATOR_FALSE_EASTING = 500000.0;
const double MERCATOR_FALSE_NORTHING = -5300000.0;
const double MERCATOR_SCALE_FACTOR = 0.9993;
const double CENTRAL_MERIDIAN = 19.0;
const double CENTRAL_PARALLEL = 52.0;

 ///////////////////////////////////////////////////////////////////////
// parallel MPI environment ///////////////////////////////////////////

int rankMPI, sizeMPI;

 ///////////////////////////////////////////////////////////////////////
// helper routines and wrappers ///////////////////////////////////////

// C++ style wrapper for Shapefile
class ShapeFileHandle {
	SHPHandle handle;
public:
	ShapeFileHandle(const char* filepath) {
		handle = SHPOpen(filepath, "rb");
		if (!handle) {
			throw std::runtime_error("file does not exist");
		}
	}
	~ShapeFileHandle(void) {
		if (handle) SHPClose(handle);
	}
	ShapeObjectHandle readObject(int index) {
		return ShapeObjectHandle(
			SHPReadObject(handle, index),
			SHPDestroyObject
		);
	}
};

// routine to read vertices from Shapefile to boost-style polygon
std::list<Polygon> readShapefile(const char* filepath) {
	ShapeFileHandle handle(filepath);

	ShapeObjectHandle shape = handle.readObject(0);
	if (!shape || shape->nSHPType != SHPT_POLYGON) {
		throw std::runtime_error("invalid file contents");
	}

	std::list<Polygon> polygons;
	const double *x = shape->padfX, *y = shape->padfY;
	for (int part=0; part<shape->nParts; ++part) {
		int vStart = (shape->nParts > 1) ? shape->panPartStart[part] : 0;
		int vEnd = (shape->nParts > 1 && part+1<shape->nParts) ? shape->panPartStart[part+1] : shape->nVertices;

		Polygon polygon;
		for (int v=vStart; v<vEnd; ++v) {
			boost::geometry::append(polygon, Point(x[v]-MERCATOR_FALSE_EASTING, y[v]-MERCATOR_FALSE_NORTHING));
		}
		polygons.push_back(std::move(polygon));
	}
	return polygons;
}

// wrapper for projection classes
class Projection {
public:
	virtual void Forward(double lat, double lon, double& x, double& y) const =0;
	virtual void Reverse(double x, double y, double& lat, double& lon) const =0;
	
	void addPointToPolygonArea(PolygonArea& polygonArea, double x, double y) const {
		double lat, lon;
		Reverse(x, y, lat, lon);
		polygonArea.AddPoint(lat, lon);
	}
};

// wrapper for transverse mercator projection
class TransMercProjection : public Projection {
	GeographicLib::TransverseMercator internal;
public:
	TransMercProjection(const Geodesic& geodesic)
	: internal(geodesic.MajorRadius(), geodesic.Flattening(), MERCATOR_SCALE_FACTOR) { }

	void Forward(double lat, double lon, double& x, double& y) const {
		internal.Forward(CENTRAL_MERIDIAN, lat, lon, x, y);
	}
	void Reverse(double x, double y, double& lat, double& lon) const {
		internal.Reverse(CENTRAL_MERIDIAN, x, y, lat, lon);
	}
};

// wrapper for Gnomonic projection
class GnomonicProjection : public Projection {
	GeographicLib::Gnomonic internal;
public:
	GnomonicProjection(const Geodesic& geodesic)
	: internal(geodesic) { }

	void Forward(double lat, double lon, double& x, double& y) const {
		internal.Forward(CENTRAL_PARALLEL, CENTRAL_MERIDIAN, lat, lon, x, y);
	}
	void Reverse(double x, double y, double& lat, double& lon) const {
		internal.Reverse(CENTRAL_PARALLEL, CENTRAL_MERIDIAN, x, y, lat, lon);
	}
};

// wrapper for Albers equal area projection
class EqualAreaProjection : public Projection {
	GeographicLib::AlbersEqualArea internal;
public:
	EqualAreaProjection(const Geodesic& geodesic)
	: internal(geodesic.MajorRadius(), geodesic.Flattening(), CENTRAL_PARALLEL, 1.0) { }

	void Forward(double lat, double lon, double& x, double& y) const {
		internal.Forward(CENTRAL_MERIDIAN, lat, lon, x, y);
	}
	void Reverse(double x, double y, double& lat, double& lon) const {
		internal.Reverse(CENTRAL_MERIDIAN, x, y, lat, lon);
	}
};

// translator from one projection to another
Polygon translate(const Projection& srcProjection, const Projection& dstProjection, const Polygon& srcPolygon) {
	Polygon dstPolygon;
	boost::geometry::for_each_point(srcPolygon, [&srcProjection, &dstProjection, &dstPolygon](const Point& point) {
		double x, y, lat, lon;
		srcProjection.Reverse(point.x(), point.y(), lat, lon);
		dstProjection.Forward(lat, lon, x, y);
		boost::geometry::append(dstPolygon, Point(x, y));
	});
	return dstPolygon;
}

 ///////////////////////////////////////////////////////////////////////
// main routines for computation of the centre ////////////////////////
void compute(const std::list<Polygon>& polygons, const Geocentric& geocentric, const Projection& projection, int N) {
	double xMin = INFINITY, xMax = -INFINITY;
	double yMin = INFINITY, yMax = -INFINITY;
	for (const Polygon& polygon : polygons) {
		Rectangle bounds;
		boost::geometry::envelope(polygon, bounds);
		xMin = std::min(xMin, bounds.min_corner().get<0>());
		xMax = std::max(xMax, bounds.max_corner().get<0>());
		yMin = std::min(yMin, bounds.min_corner().get<1>());
		yMax = std::max(yMax, bounds.max_corner().get<1>());
	}

	double x[N+1], y[N+1];
	for (int i=0; i<=N; ++i) {
		x[i] = xMin + (xMax - xMin) * i/N;
		y[i] = yMin + (yMax - yMin) * i/N;
	}

	double X0, Y0, Z0;
	geocentric.Forward(CENTRAL_PARALLEL, CENTRAL_MERIDIAN, 0.0, X0, Y0, Z0);

	Polygon cell;
	std::list<Polygon> intersection;

	double areaNode = 0;
	double XNode = 0, YNode = 0, ZNode = 0;

	MPI_Barrier(MPI_COMM_WORLD);
	for (int ix=rankMPI; ix<N; ix+=sizeMPI) {
		for (int iy=0; iy<N; ++iy) {
			cell.clear();
			boost::geometry::append(cell, Point(x[ix], y[iy]));
			boost::geometry::append(cell, Point(x[ix], y[iy+1]));
			boost::geometry::append(cell, Point(x[ix+1], y[iy+1]));
			boost::geometry::append(cell, Point(x[ix+1], y[iy]));
			boost::geometry::append(cell, Point(x[ix], y[iy]));

			for (const Polygon& polygon : polygons) {
				intersection.clear();
				boost::geometry::intersection(polygon, cell, intersection);
				for (const Polygon& part : intersection) {
					Point centroid(0, 0);
					double lat, lon, X, Y, Z;
					boost::geometry::centroid(part, centroid);
					projection.Reverse(centroid.x(), centroid.y(), lat, lon);
					geocentric.Forward(lat, lon, 0.0, X, Y, Z);

					double area = boost::geometry::area(part);
					XNode += (X - X0) * area;
					YNode += (Y - Y0) * area;
					ZNode += (Z - Z0) * area;

					// summation of the area
					areaNode += area;
				}
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	double X, Y, Z, area;
	MPI_Reduce(&areaNode, &area, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&XNode, &X, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&YNode, &Y, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&ZNode, &Z, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (!rankMPI) {
		X = X0 + X / area;
		Y = Y0 + Y / area;
		Z = Z0 + Z / area;
		double lat, lon, h;
		geocentric.Reverse(X, Y, Z, lat, lon, h);
		double latD, latM, latS, lonD, lonM, lonS;
		DMS::Encode(lat, latD, latM, latS);
		DMS::Encode(lon, lonD, lonM, lonS);
		printf("  %10.3f %12.10f %12.10f %.3f %.0f°%02.0f′%05.2f″ %.0f°%02.0f′%05.2f″",
			1.0e-6*area, lat, lon, h, latD, latM, latS, lonD, lonM, lonS);
		fflush(stdout);
	}
}

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rankMPI);
	MPI_Comm_size(MPI_COMM_WORLD, &sizeMPI);
	if (argc < 2) {
		if (!rankMPI) {
			fprintf(stderr, "USAGE: %s N\n", argv[0]);
		}
		return 1;
	}
	try {
		const int N = atoi(argv[1]);
		if (N < 10) {
			throw std::runtime_error("N must be at least 10");
		}

		Geodesic geodesic(GRS80_A, GRS80_F);
		Geocentric geocentric(GRS80_A, GRS80_F);
		TransMercProjection transMerc(geodesic);
		EqualAreaProjection equalArea(geodesic);
		if (!rankMPI) {
			printf("%6d", N);
			fflush(stdout);
		}

		std::list<Polygon> polygonsInTransMerc = readShapefile("data/Państwo.shp");
		std::list<Polygon> polygonsInEqualArea;
		for (const Polygon& polygon : polygonsInTransMerc) {
			polygonsInEqualArea.push_back(translate(transMerc, equalArea, polygon));
		}
		compute(polygonsInEqualArea, geocentric, equalArea, N);
		if (!rankMPI) {
			putchar('\n');
		}

	} catch (const std::exception& ex) {
		if (!rankMPI) {
			fprintf(stderr, "ERROR: %s\n", ex.what());
		}
		return 1;
	}
	MPI_Finalize();
}
