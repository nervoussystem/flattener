/*
  Copyright (c) 2009 Erin Catto http://www.box2d.org
  Copyright (c) 2016-2017 Lester Hedges <lester.hedges+aabbcc@gmail.com>

  This software is provided 'as-is', without any express or implied
  warranty. In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.

  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.

  3. This notice may not be removed or altered from any source distribution.

  This code was adapted from parts of the Box2D Physics Engine,
  http://www.box2d.org
*/

#ifndef _AABB_H
#define _AABB_H

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

/// Null node flag.
const unsigned int NULL_NODE = 0xffffffff;

namespace aabb
{
    /*! \brief The axis-aligned bounding box object.

        Axis-aligned bounding boxes (AABBs) store information for the minimum
        orthorhombic bounding-box for an object in two- or three-dimensional
        space (the bounding box is either a rectangle, or rectangular prism).

        Class member functions provide functionality for merging AABB objects
        and testing overlap with other AABBs.
     */
    class AABB
    {
    public:
        /// Constructor.
        AABB();

        //! Constructor.
        /*! \param dimension
                The dimensionality of the system.
         */
        AABB(unsigned int);

        //! Constructor.
        /*! \param lowerBound_
                The lower bound in each dimension.

            \param upperBound_
                The upper bound in each dimension.
         */
        AABB(const std::vector<double>&, const std::vector<double>&);

        /// Compute the surface area of the box.
        double computeSurfaceArea() const;

        /// Get the surface area of the box.
        double getSurfaceArea() const;

        //! Merge two AABBs into this one.
        /*! \param aabb1
                A reference to the first AABB.

            \param aabb2
                A reference to the second AABB.
         */
        void merge(const AABB&, const AABB&);

        //! Test whether the AABB is contained within this one.
        /*! \param aabb
                A reference to the AABB.

            \return
                Whether the AABB is fully contained.
         */
        bool contains(const AABB&) const;

        //! Test whether the AABB overlaps this one.
        /*! \param aabb
                A reference to the AABB.

            \return
                Whether the AABB overlaps.
         */
        bool overlaps(const AABB&) const;

        //! Compute the centre of the AABB.
        /*! \returns
                The position vector of the AABB centre.
         */
        std::vector<double> computeCentre();

        //! Set the dimensionality of the AABB.
        /*! \param dimension
                The dimensionality of the system.
         */
        void setDimension(unsigned int);

        /// Lower bound of AABB in each dimension.
        std::vector<double> lowerBound;

        /// Upper bound of AABB in each dimension.
        std::vector<double> upperBound;

        /// The position of the AABB centre.
        std::vector<double> centre;

        /// The AABB's surface area.
        double surfaceArea;
    };

    /*! \brief A node of the AABB tree.

        Each node of the tree contains an AABB object which corresponds to a
        particle, or a group of particles, in the simulation box. The AABB
        objects of individual particles are "fattened" before they are stored
        to avoid having to continually update and rebalance the tree when
        displacements are small.

        Nodes are aware of their position within in the tree. The isLeaf member
        function allows the tree to query whether the node is a leaf, i.e. to
        determine whether it holds a single particle.
     */
    struct Node
    {
        /// Constructor.
        Node();

        /// The fattened axis-aligned bounding box.
        AABB aabb;

		int index;

		std::vector<Node *> children;
        /// Index of the parent node.
        unsigned int parent;

        bool isLeaf() const;
    };

    /*! \brief The dynamic AABB tree.

        The dynamic AABB tree is a hierarchical data structure that can be used
        to efficiently query overlaps between objects of arbitrary shape and
        size that lie inside of a simulation box. Support is provided for
        periodic and non-periodic boxes, as well as boxes with partial
        periodicity, e.g. periodic along specific axes.
     */
    class Tree
    {
    public:
		Tree();

		void build(std::vector<AABB> & boxes);
		std::vector<unsigned int> Tree::query(const AABB& aabb);

		void splitTree(Node * node);
		//template<typename T>
		//std::pair<std::array<double, 3>, int> closestPointAndPrimitive(const std::array<double, 3> & pt, double estimatedRadius, T & closestOperator);

		Node * root;

		template<typename T>
		void closestPointAndPrimitiveH(AABB & aabb, const std::array<double, 3> & pt, Node * n, double &closestD, std::array<double, 3> & closestPt, int & closestI, T & closestOperator) {
			std::array<double, 3> tPt;
			for (int i = 0; i < n->children.size(); ++i) {
				Node * child = n->children[i];
				if (aabb.overlaps(child->aabb)) {
					if (child->index >= 0) {
						double d = closestOperator(pt, child->index, tPt);
						if (d < closestD) {
							closestD = d;
							float dsq = sqrt(d);
							closestI = child->index;
							closestPt = tPt;
							if (d == 0) {
								return;
							}
							else {
								for (int j = 0; j < 3; ++j) {
									aabb.lowerBound[j] = pt[j] - dsq;
									aabb.upperBound[j] = pt[j] + dsq;
								}
							}
						}
					}
					else {
						closestPointAndPrimitiveH(aabb, pt,child, closestD, closestPt, closestI, closestOperator);
					}
				}
			}
		}

		template<typename T>
		std::pair<std::array<double, 3>, int> closestPointAndPrimitive(const std::array<double, 3> & pt, double estimatedRadius, T & closestOperator) {
			AABB q(3);
			for (int i = 0; i < 3; ++i) {
				q.lowerBound[i] = pt[i] - estimatedRadius;
				q.upperBound[i] = pt[i] + estimatedRadius;
			}
			double closestD = estimatedRadius*estimatedRadius;// std::numeric_limits<double>::max();
			int closestI = 0;
			std::array<double, 3> closestPt;
			closestPointAndPrimitiveH(q, pt, root, closestD, closestPt, closestI, closestOperator);
			return make_pair(closestPt, closestI);
		}
    };
}

#endif /* _AABB_H */
