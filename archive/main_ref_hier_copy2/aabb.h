#ifndef AABB_H
#define AABB_H

#include "../Eigen/Dense"
#include <limits>

#include "ray.h"

const double infinity = std::numeric_limits<double>::infinity();

class interval {
  public:
    double min, max;

    interval() : min(+infinity), max(-infinity) {} // Default interval is empty

    interval(double min, double max) : min(min), max(max) {}

    interval(const interval& a, const interval& b) {
        // Create the interval tightly enclosing the two input intervals.
        min = a.min <= b.min ? a.min : b.min;
        max = a.max >= b.max ? a.max : b.max;
    }

    double size() const {
        return max - min;
    }

    bool contains(double x) const {
        return min <= x && x <= max;
    }

    bool surrounds(double x) const {
        return min < x && x < max;
    }

    double clamp(double x) const {
        if (x < min) return min;
        if (x > max) return max;
        return x;
    }

    interval expand(double delta) const {
        auto padding = delta/2;
        return interval(min - padding, max + padding);
    }

    static const interval empty, universe;
};

const interval interval::empty    = interval(+infinity, -infinity);
const interval interval::universe = interval(-infinity, +infinity);

interval operator+(const interval& ival, double displacement) {
    return interval(ival.min + displacement, ival.max + displacement);
}

interval operator+(double displacement, const interval& ival) {
    return ival + displacement;
}

// Intersect function check if the specific ray intersects the box,
// We don't need it to return the information about the intersection,
// as we don't care about it.
class aabb {
  public:
    interval x, y, z;

    aabb() {} // The default AABB is empty, since intervals are empty by default.

    aabb(const interval& x, const interval& y, const interval& z)
      : x(x), y(y), z(z)
    {
        pad_to_minimums();
    }

    aabb(const Eigen::Vector3d &a, const Eigen::Vector3d& b) {
        // Treat the two points a and b as extrema for the bounding box, so we don't require a
        // particular minimum/maximum coordinate order.

        x = (a[0] <= b[0]) ? interval(a[0], b[0]) : interval(b[0], a[0]);
        y = (a[1] <= b[1]) ? interval(a[1], b[1]) : interval(b[1], a[1]);
        z = (a[2] <= b[2]) ? interval(a[2], b[2]) : interval(b[2], a[2]);

        pad_to_minimums();
    }

    aabb(const aabb& box0, const aabb& box1) {
        x = interval(box0.x, box1.x);
        y = interval(box0.y, box1.y);
        z = interval(box0.z, box1.z);
    }

    const interval& axis_interval(int n) const {
        if (n == 1) return y;
        if (n == 2) return z;
        return x;
    }

    bool intersect(const Ray &r, interval ray_t) const {
        const Eigen::Vector3d ray_orig = r.o;
        const Eigen::Vector3d   ray_dir  = r.d;

        for (int axis = 0; axis < 3; axis++) {
            const interval& ax = axis_interval(axis);
            const double adinv = 1.0 / ray_dir[axis];

            auto t0 = (ax.min - ray_orig[axis]) * adinv;
            auto t1 = (ax.max - ray_orig[axis]) * adinv;

            if (t0 < t1) {
                if (t0 > ray_t.min) ray_t.min = t0;
                if (t1 < ray_t.max) ray_t.max = t1;
            } else {
                if (t1 > ray_t.min) ray_t.min = t1;
                if (t0 < ray_t.max) ray_t.max = t0;
            }

            if (ray_t.max <= ray_t.min)
                return false;
        }
        return true;
    }

    int longest_axis() const {
        // Returns the index of the longest axis of the bounding box.

        if (x.size() > y.size())
            return x.size() > z.size() ? 0 : 2;
        else
            return y.size() > z.size() ? 1 : 2;
    }

    static const aabb empty, universe;

  private:

    void pad_to_minimums() {
        // Adjust the AABB so that no side is narrower than some delta, padding if necessary.

        double delta = 0.0001;
        if (x.size() < delta) x = x.expand(delta);
        if (y.size() < delta) y = y.expand(delta);
        if (z.size() < delta) z = z.expand(delta);
    }
};

const aabb aabb::empty    = aabb(interval::empty,    interval::empty,    interval::empty);
const aabb aabb::universe = aabb(interval::universe, interval::universe, interval::universe);

aabb operator+(const aabb& bbox, const Eigen::Vector3d &offset) {
    return aabb(bbox.x + offset.x(), bbox.y + offset.y(), bbox.z + offset.z());
}

aabb operator+(const Eigen::Vector3d &offset, const aabb& bbox) {
    return bbox + offset;
}


// class bvh_node : public intersectable {
//   public:
//     bvh_node(intersectable_list list) : bvh_node(list.objects, 0, list.size()) {
//         // There's a C++ subtlety here. This constructor (without span indices) creates an
//         // implicit copy of the hittable list, which we will modify. The lifetime of the copied
//         // list only extends until this constructor exits. That's OK, because we only need to
//         // persist the resulting bounding volume hierarchy.
//     }

//     bvh_node(const std::vector<std::shared_ptr<intersectable>>& objects, size_t start, size_t end) {
//         // Build the bounding box of the span of source objects.
//         bbox = aabb::empty;
//         for (size_t object_index=start; object_index < end; object_index++)
//             bbox = aabb(bbox, objects[object_index]->bbox);

//         int axis = bbox.longest_axis();

//         auto comparator = (axis == 0) ? box_x_compare
//                         : (axis == 1) ? box_y_compare
//                                       : box_z_compare;

//         size_t object_span = end - start;

//         if (object_span == 1) {
//             left = right = objects[start];
//         } else if (object_span == 2) {
//             left = objects[start];
//             right = objects[start+1];
//         } else {
//             std::sort(std::begin(objects) + start, std::begin(objects) + end, comparator);

//             auto mid = start + object_span/2;
//             left = std::make_shared<bvh_node>(objects, start, mid);
//             right = std::make_shared<bvh_node>(objects, mid, end);
//         }
//     }

//   private:
//     std::shared_ptr<intersectable> left;
//     std::shared_ptr<intersectable> right;
//     aabb bbox;

//     static bool box_compare(
//         const std::shared_ptr<intersectable> a, const std::shared_ptr<intersectable> b, int axis_index
//     ) {
//         auto a_axis_interval = a->bbox.axis_interval(axis_index);
//         auto b_axis_interval = b->bbox.axis_interval(axis_index);
//         return a_axis_interval.min < b_axis_interval.min;
//     }

//     static bool box_x_compare (const std::shared_ptr<intersectable> a, const std::shared_ptr<intersectable> b) {
//         return box_compare(a, b, 0);
//     }

//     static bool box_y_compare (const std::shared_ptr<intersectable> a, const std::shared_ptr<intersectable> b) {
//         return box_compare(a, b, 1);
//     }

//     static bool box_z_compare (const std::shared_ptr<intersectable> a, const std::shared_ptr<intersectable> b) {
//         return box_compare(a, b, 2);
//     }

    
// };

#endif
