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
typedef boost::geometry::model::ring<Point> Ring;
typedef boost::geometry::model::polygon<Point> Polygon;
typedef boost::geometry::model::multi_polygon<Polygon> Shape;
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
Shape readShapefile(const char* filepath, int index) {
	ShapeFileHandle handle(filepath);

	ShapeObjectHandle shape = handle.readObject(index);
	if (!shape || shape->nSHPType != SHPT_POLYGON) {
		throw std::runtime_error("invalid file contents");
	}

	Shape result;
	std::list<Ring> innerRings;

	const double *x = shape->padfX, *y = shape->padfY;
	for (int part=0; part<shape->nParts; ++part) {
		int vStart = (shape->nParts > 1) ? shape->panPartStart[part] : 0;
		int vEnd = (shape->nParts > 1 && part+1<shape->nParts) ? shape->panPartStart[part+1] : shape->nVertices;

		Ring ring;
		for (int v=vStart; v<vEnd; ++v) {
			boost::geometry::append(ring, Point(x[v]-MERCATOR_FALSE_EASTING, y[v]-MERCATOR_FALSE_NORTHING));
		}
		if (boost::geometry::area(ring) < 0) {
			innerRings.push_back(ring);
		} else {
			result.push_back({ring});
		}
	}

	for (const Ring& ring : innerRings) {
		bool found = false;
		for (Polygon& polygon : result) {
			if (boost::geometry::within(ring, polygon)) {
				if (found) {
					throw std::runtime_error("found inner ring with multiple matching outer rings");
				}
				polygon.inners().push_back(ring);
				found = true;
			}
		}
		if (!found) {
			throw std::runtime_error("found inner ring with no matching outer ring");
		}
	}

	return result;
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

 ///////////////////////////////////////////////////////////////////////
// translators from one projection to another /////////////////////////
void translateRing(const Projection& srcProjection, const Projection& dstProjection, Ring& ring) {
	boost::geometry::for_each_point(ring, [&srcProjection, &dstProjection](Point& point) {
		double x, y, lat, lon;
		srcProjection.Reverse(point.x(), point.y(), lat, lon);
		dstProjection.Forward(lat, lon, x, y);
		point.x(x); point.y(y);
	});
}

void translatePolygon(const Projection& srcProjection, const Projection& dstProjection, Polygon& polygon) {
	translateRing(srcProjection, dstProjection, polygon.outer());
	for (auto& inner : polygon.inners()) {
		translateRing(srcProjection, dstProjection, inner);
	}
}

void translateShape(const Projection& srcProjection, const Projection& dstProjection, Shape& shape) {
	for (auto& polygon : shape) {
		translatePolygon(srcProjection, dstProjection, polygon);
	}
}

 ///////////////////////////////////////////////////////////////////////
// main routines for computation of the centre ////////////////////////
void compute(const Shape& shape, const Geocentric& geocentric, const Projection& projection, int N) {
	Rectangle bounds;
	boost::geometry::envelope(shape, bounds);
	double xMin = bounds.min_corner().x();
	double xMax = bounds.max_corner().x();
	double yMin = bounds.min_corner().y();
	double yMax = bounds.max_corner().y();

	double x[N+1], y[N+1];
	for (int i=0; i<=N; ++i) {
		x[i] = xMin + (xMax - xMin) * i/N;
		y[i] = yMin + (yMax - yMin) * i/N;
	}

	double X0, Y0, Z0;
	geocentric.Forward(CENTRAL_PARALLEL, CENTRAL_MERIDIAN, 0.0, X0, Y0, Z0);

	Ring cell;
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

			intersection.clear();
			boost::geometry::intersection(shape, cell, intersection);
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
		printf(" %10.3f %12.10f %12.10f %.3f %.0f°%02.0f′%05.2f″ %.0f°%02.0f′%05.2f″",
			1.0e-6*area, lat, lon, h, latD, latM, latS, lonD, lonM, lonS);
		fflush(stdout);
	}
}

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rankMPI);
	MPI_Comm_size(MPI_COMM_WORLD, &sizeMPI);
	if (argc < 3) {
		if (!rankMPI) {
			fprintf(stderr, "USAGE: %s N path_to_shapefile [ object_index ]\n", argv[0]);
		}
		return 1;
	}
	try {
		const int N = atoi(argv[1]);
		const char* path = argv[2];
		const int index = argc>3 ? atoi(argv[3]) : 0;
		if (N < 10) {
			throw std::runtime_error("N must be at least 10");
		}

		Geodesic geodesic(GRS80_A, GRS80_F);
		Geocentric geocentric(GRS80_A, GRS80_F);
		TransMercProjection transMerc(geodesic);
		EqualAreaProjection equalArea(geodesic);
		if (!rankMPI) {
			printf("%6d ", N);
			fflush(stdout);
		}

		Shape shape = readShapefile(path, index);
		translateShape(transMerc, equalArea, shape);
		compute(shape, geocentric, equalArea, N);
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
