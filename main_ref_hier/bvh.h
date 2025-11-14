#ifndef BVH_H
#define BVH_H

#include "../Eigen/Dense"
#include <limits>

#include "ray.h"
#include "intersectable.h"
#include "aabb.h"

class bvh_node : public intersectable {
  public:
    bvh_node(intersectable_list list) : bvh_node(list.objects, 0, list.objects.size()) {
        // There's a C++ subtlety here. This constructor (without span indices) creates an
        // implicit copy of the hittable list, which we will modify. The lifetime of the copied
        // list only extends until this constructor exits. That's OK, because we only need to
        // persist the resulting bounding volume hierarchy.
    }

    bvh_node(std::vector<std::shared_ptr<intersectable>>& objects, size_t start, size_t end) {
        // Build the bounding box of the span of source objects.
        bbox = aabb::empty;
        for (size_t object_index=start; object_index < end; object_index++)
            bbox = aabb(bbox, objects[object_index]->bounding_box());
        // std::cout << " Object size: "<< start << " and " <<  end << " and size " << objects.size() << ". \n";

        int axis = bbox.longest_axis();

        auto comparator = (axis == 0) ? box_x_compare
                        : (axis == 1) ? box_y_compare
                                      : box_z_compare;

        size_t object_span = end - start;

        if (object_span == 1) {
            left = right = objects[start];
        } else if (object_span == 2) {
            left = objects[start];
            right = objects[start+1];
        } else {
            std::sort(std::begin(objects) + start, std::begin(objects) + end, comparator);

            auto mid = start + object_span/2;
            left = std::make_shared<bvh_node>(objects, start, mid);
            right = std::make_shared<bvh_node>(objects, mid, end);
        }
    }

    bool intersect(const Ray& r, interval ray_t, intersect_record& rec) const override {
        if (!bbox.intersect(r, ray_t))
            return false;

        bool hit_left = left->intersect(r, ray_t, rec);
        bool hit_right = right->intersect(r, interval(ray_t.min, hit_left ? rec.t : ray_t.max), rec);

        return hit_left || hit_right;
    }


  private:
    std::shared_ptr<intersectable> left;
    std::shared_ptr<intersectable> right;

    static bool box_compare(
        const std::shared_ptr<intersectable> a, const std::shared_ptr<intersectable> b, int axis_index
    ) {
        auto a_axis_interval = a->bounding_box().axis_interval(axis_index);
        auto b_axis_interval = b->bounding_box().axis_interval(axis_index);
        return a_axis_interval.min < b_axis_interval.min;
    }

    static bool box_x_compare (const std::shared_ptr<intersectable> a, const std::shared_ptr<intersectable> b) {
        return box_compare(a, b, 0);
    }

    static bool box_y_compare (const std::shared_ptr<intersectable> a, const std::shared_ptr<intersectable> b) {
        return box_compare(a, b, 1);
    }

    static bool box_z_compare (const std::shared_ptr<intersectable> a, const std::shared_ptr<intersectable> b) {
        return box_compare(a, b, 2);
    }
};

#endif
