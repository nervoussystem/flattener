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

#include "AABB1.h"

namespace aabb
{
	AABB::AABB()
	{
	}

	AABB::AABB(unsigned int dimension)
	{
		assert((dimension == 2) || (dimension == 3));

		lowerBound.resize(dimension);
		upperBound.resize(dimension);

		for (int i = 0; i < dimension; ++i) {
			lowerBound[i] = std::numeric_limits<double>::max();
			upperBound[i] = std::numeric_limits<double>::min();
		}
	}

	AABB::AABB(const std::vector<double>& lowerBound_, const std::vector<double>& upperBound_) :
		lowerBound(lowerBound_), upperBound(upperBound_)
	{
		surfaceArea = computeSurfaceArea();
		centre = computeCentre();
	}

	double AABB::computeSurfaceArea() const
	{
		// Calculate the perimeter of the 2D AABB.
		if (lowerBound.size() == 2)
		{
			double wx = upperBound[0] - lowerBound[0];
			double wy = upperBound[1] - lowerBound[1];
			return 2.0 * (wx + wy);
		}

		// Calculate the surface area of the 3D AABB.
		else
		{
			double wx = upperBound[0] - lowerBound[0];
			double wy = upperBound[1] - lowerBound[1];
			double wz = upperBound[2] - lowerBound[2];
			return 2.0 * (wx*wy + wx*wz + wy*wz);
		}
	}

	double AABB::getSurfaceArea() const
	{
		return surfaceArea;
	}

	void AABB::merge(const AABB& aabb1, const AABB& aabb2)
	{
		assert(aabb1.lowerBound.size() == aabb2.lowerBound.size());
		assert(aabb1.upperBound.size() == aabb2.upperBound.size());

		lowerBound.resize(aabb1.lowerBound.size());
		upperBound.resize(aabb1.lowerBound.size());

		for (unsigned int i = 0; i < lowerBound.size(); i++)
		{
			lowerBound[i] = std::min(aabb1.lowerBound[i], aabb2.lowerBound[i]);
			upperBound[i] = std::max(aabb1.upperBound[i], aabb2.upperBound[i]);
		}

		surfaceArea = computeSurfaceArea();
		centre = computeCentre();
	}

	bool AABB::contains(const AABB& aabb) const
	{
		assert(aabb.lowerBound.size() == lowerBound.size());

		for (unsigned int i = 0; i < lowerBound.size(); i++)
		{
			if (aabb.lowerBound[i] < lowerBound[i]) return false;
			if (aabb.upperBound[i] > upperBound[i]) return false;
		}

		return true;
	}

	bool AABB::overlaps(const AABB& aabb) const
	{
		assert(aabb.lowerBound.size() == lowerBound.size());

		if (lowerBound.size() == 2)
		{
			return !(aabb.upperBound[0] < lowerBound[0]
				|| aabb.lowerBound[0] > upperBound[0]
				|| aabb.upperBound[1] < lowerBound[1]
				|| aabb.lowerBound[1] > upperBound[1]
				);
		}
		else
		{
			return !(aabb.upperBound[0] < lowerBound[0]
				|| aabb.lowerBound[0] > upperBound[0]
				|| aabb.upperBound[1] < lowerBound[1]
				|| aabb.lowerBound[1] > upperBound[1]
				|| aabb.upperBound[2] < lowerBound[2]
				|| aabb.lowerBound[2] > upperBound[2]
				);
		}
	}

	std::vector<double> AABB::computeCentre()
	{
		std::vector<double> position(lowerBound.size());

		for (unsigned int i = 0; i < position.size(); i++)
			position[i] = 0.5 * (lowerBound[i] + upperBound[i]);

		return position;
	}

	void AABB::setDimension(unsigned int dimension)
	{
		assert((dimension == 2) || (dimension == 3));

		lowerBound.resize(dimension);
		upperBound.resize(dimension);
	}

	Node::Node() :
		aabb(3)
	{
		index = -1;
	}

	bool Node::isLeaf() const
	{
		return index < 0;
	}

	Tree::Tree()
	{
		root = NULL;
	}

	void Tree::build(std::vector<AABB> & boxes) {
		root = new Node();
		root->aabb = boxes[0];
		for (int i = 0; i < boxes.size(); ++i) {
			Node * leaf = new Node();
			leaf->index = i;
			leaf->aabb = boxes[i];
			root->children.push_back(leaf);
			root->aabb.merge(root->aabb, leaf->aabb);
		}
		splitTree(root);
	}

	void Tree::splitTree(Node * node) {
		int  axis = 0;
		double len = node->aabb.upperBound[0] - node->aabb.lowerBound[0];
		if (node->aabb.upperBound[1] - node->aabb.lowerBound[1] > len) {
			len = node->aabb.upperBound[1] - node->aabb.lowerBound[1];
			axis = 1;
		}
		if (node->aabb.upperBound[2] - node->aabb.lowerBound[2] > len) {
			len = node->aabb.upperBound[2] - node->aabb.lowerBound[2];
			axis = 2;
		}
		sort(node->children.begin(), node->children.end(), [axis](Node * n1, Node * n2) {return n1->aabb.centre[axis] < n2->aabb.centre[axis]; });
		double mid = node->aabb.centre[axis];
		Node * leaf1 = new Node();
		Node * leaf2 = new Node();
		/*
		for (int i = 0; i < node->children.size(); ++i) {
			Node * child = node->children[i];
			if (child->aabb.centre[axis] > mid) {
				leaf1->children.push_back(child);
				leaf1->aabb.merge(leaf1->aabb, child->aabb);
			}
			else {
				leaf2->children.push_back(child);
				leaf2->aabb.merge(leaf2->aabb, child->aabb);
			}
		}
		*/
		int median = node->children.size() / 2;
		for (int i = 0; i < median; ++i) {
			Node * child = node->children[i];
			leaf1->children.push_back(child);
			leaf1->aabb.merge(leaf1->aabb, child->aabb);
		}
		for (int i = median,sz = node->children.size(); i < sz; ++i) {
			Node * child = node->children[i];
			leaf2->children.push_back(child);
			leaf2->aabb.merge(leaf1->aabb, child->aabb);
		}
		node->children.clear();
		/*
		if (leaf1->children.size() == 0) {
			node->children.insert(node->children.end(), leaf2->children.begin(), leaf2->children.end());
		}
		else if (leaf2->children.size() == 0) {
			node->children.insert(node->children.end(), leaf1->children.begin(), leaf1->children.end());
		}
		else {
		*/
			node->children.push_back(leaf1);
			node->children.push_back(leaf2);
			if (leaf1->children.size() > 3) {
				splitTree(leaf1);
			}
			if (leaf2->children.size() > 3) {
				splitTree(leaf2);
			}
		//}
	}

	void queryH(const AABB & aabb, Node * n, std::vector<unsigned int> & out) {
		for (int i = 0; i < n->children.size(); ++i) {
			Node * child = n->children[i];
			if (aabb.overlaps(child->aabb)) {
				if (child->index >= 0) {
					out.push_back(child->index);
				}
				else {
					queryH(aabb, child, out);
				}
			}
		}
	}

	std::vector<unsigned int> Tree::query(const AABB& aabb)
	{
		// Make sure the tree isn't empty.
		if (root == NULL)
		{
			return std::vector<unsigned int>();
		}
		std::vector<unsigned int> out;
		queryH(aabb, root, out);
		return out;
	}

	
}