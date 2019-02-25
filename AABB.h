#pragma once

#include <boost/iterator.hpp>
#include "ofMain.h"
#include "nsHEMesh.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>

typedef CGAL::Simple_cartesian<double> BK;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel BK;

// The triangles are stored in a flat array of indices 
// referring to an array of points: three consecutive
// indices represent a triangle.
typedef std::vector<ofIndexType>::const_iterator Point_index_iterator;
// Let us now define the iterator on triangles that the tree needs:
class Triangle_iterator
    : public boost::iterator_adaptor<
    Triangle_iterator               // Derived
    , Point_index_iterator            // Base
    , boost::use_default              // Value
    , boost::forward_traversal_tag    // CategoryOrTraversal
    >
{
public:
    Triangle_iterator()
        : Triangle_iterator::iterator_adaptor_() {}
    explicit Triangle_iterator(Point_index_iterator p,vector<ofVec3f> * cont)
        : Triangle_iterator::iterator_adaptor_(p) {point_container = cont;}
	vector<ofVec3f>* point_container;
private:

    friend class boost::iterator_core_access;
    void increment() { this->base_reference() += 3; }
};
// The following primitive provides the conversion facilities between
// my own triangle and point types and the CGAL ones
struct My_triangle_primitive {
public:
    typedef Triangle_iterator    Id;
    // the CGAL types returned
    typedef BK::Point_3    Point;
    typedef BK::Triangle_3 Datum;
    // a static pointer to the vector containing the points
    // is needed to build the triangles on the fly:
    vector<ofVec3f>* point_container;
private:
    Id m_it; // this is what the AABB tree stores internally
public:
    My_triangle_primitive() {} // default constructor needed
    // the following constructor is the one that receives the iterators from the 
    // iterator range given as input to the AABB_tree
    My_triangle_primitive(Triangle_iterator a)
		: m_it(a) {point_container = a.point_container;}
    Id id() const { return m_it; }
    // on the fly conversion from the internal data to the CGAL types
    Datum datum() const
    { 
		Point_index_iterator p_it = m_it.base();
        const ofVec3f& mp = (*point_container)[*p_it];
        Point p(mp.x, mp.y, mp.z);
        ++p_it;
        const ofVec3f& mq = (*point_container)[*p_it];
        Point q(mq.x, mq.y, mq.z);
        ++p_it;
        const ofVec3f& mr = (*point_container)[*p_it];
        Point r(mr.x, mr.y, mr.z);
        return Datum(p, q, r); // assembles triangle from three points
    }
    // one point which must be on the primitive
    Point reference_point() const
    { 
        const ofVec3f& mp = (*point_container)[*m_it];
        return Point(mp.x, mp.y, mp.z);
    }
};

typedef vector<nsHEFace *>::iterator HEFace_iterator;

struct HEFace_primitive {
public:
	typedef HEFace_iterator    Id;
	// the CGAL types returned
	typedef BK::Point_3    Point;
	typedef BK::Triangle_3 Datum;
	// a static pointer to the vector containing the points
	// is needed to build the triangles on the fly:
private:
	Id m_it; // this is what the AABB tree stores internally
public:
	HEFace_primitive() {} // default constructor needed
							   // the following constructor is the one that receives the iterators from the 
							   // iterator range given as input to the AABB_tree
	HEFace_primitive(HEFace_iterator a)
		: m_it(a) {
	}
	Id id() const { return m_it; }
	// on the fly conversion from the internal data to the CGAL types
	Datum datum() const
	{
		nsHEdge * e = (*m_it)->edge;
		const ofVec3f& mp = e->vertex->getPosition();
		Point p(mp.x, mp.y, mp.z);
		e = e->next;
		const ofVec3f& mq = e->vertex->getPosition();
		Point q(mq.x, mq.y, mq.z);
		e = e->next;
		const ofVec3f& mr = e->vertex->getPosition();
		Point r(mr.x, mr.y, mr.z);
		return Datum(p, q, r); // assembles triangle from three points
	}
	// one point which must be on the primitive
	Point reference_point() const
	{
		const ofVec3f& mp = (*m_it)->edge->vertex->getPosition();
		return Point(mp.x, mp.y, mp.z);
	}
};

typedef list<nsHEEdge *>::iterator HEEdge_iterator;

struct HEEdge_primitive {
public:
	typedef HEEdge_iterator    Id;
	// the CGAL types returned
	typedef BK::Point_3    Point;
	typedef BK::Segment_3 Datum;
	// a static pointer to the vector containing the points
	// is needed to build the triangles on the fly:
private:
	Id m_it; // this is what the AABB tree stores internally
public:
	HEEdge_primitive() {} // default constructor needed
						  // the following constructor is the one that receives the iterators from the 
						  // iterator range given as input to the AABB_tree
	HEEdge_primitive(HEEdge_iterator a)
		: m_it(a) {
	}
	Id id() const { return m_it; }
	// on the fly conversion from the internal data to the CGAL types
	Datum datum() const
	{
		nsHEdge * e = (*m_it)->edge;
		const ofVec3f& mp = e->vertex->getPosition();
		Point p(mp.x, mp.y, mp.z);
		e = e->pair;
		const ofVec3f& mq = e->vertex->getPosition();
		Point q(mq.x, mq.y, mq.z);
		return Datum(p, q); // assembles triangle from three points
	}
	// one point which must be on the primitive
	Point reference_point() const
	{
		const ofVec3f& mp = (*m_it)->edge->vertex->getPosition();
		return Point(mp.x, mp.y, mp.z);
	}
};

// types
typedef CGAL::AABB_traits<BK, My_triangle_primitive> My_AABB_traits;
typedef CGAL::AABB_tree<My_AABB_traits> BTree;

typedef BTree::Object_and_primitive_id BObject_and_primitive_id;
typedef BTree::Point_and_primitive_id BPoint_and_primitive_id;
typedef BTree::Primitive_id BPrimitive_id;

typedef CGAL::AABB_traits<BK, HEFace_primitive> HE_AABB_traits;
typedef CGAL::AABB_tree<HE_AABB_traits> HTree;

typedef HTree::Object_and_primitive_id HObject_and_primitive_id;
typedef HTree::Point_and_primitive_id HPoint_and_primitive_id;
typedef HTree::Primitive_id HPrimitive_id;

typedef BK::Ray_3 BRay;
typedef BK::Point_3 BPoint;
typedef BK::Segment_3 BSegment;

typedef CGAL::AABB_traits<BK, HEEdge_primitive> STraits;
typedef CGAL::AABB_tree<STraits> STree;


typedef boost::optional< STree::Intersection_and_primitive_id<BK::Triangle_3>::Type > Seg_triangle_intersection;

typedef boost::optional< HTree::Intersection_and_primitive_id<BK::Segment_3>::Type > Segment_intersection;
typedef boost::optional< HTree::Intersection_and_primitive_id<BK::Triangle_3>::Type > Triangle_intersection;
