//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    Intersection inter = intersect(ray);

    if(!inter.happened)
        return Vector3f(0, 0, 0);

    if(inter.m->hasEmission()){
        if(depth == 0)
            return inter.m->hasEmission();
        else
            return Vector3f(0, 0, 0);
    }

    Vector3f L_dir(0, 0, 0);
    Vector3f L_indir(0, 0, 0);

    //direct
    // 随机 sample 灯光，用该 sample 的结果判断射线是否击中光源
    Intersection light_sample;
    float pdf;
    sampleLight(light_sample, pdf);
    Vector3f obj_point_light = normalize(light_sample.coords - inter.coords);

    Intersection inter_light = intersect(Ray(inter.coords, obj_point_light));
    if(inter_light.happened && (inter_light.coords - light_sample.coords).norm() < 1e-2){
        Vector3f f_r = inter.m->eval(ray.direction, obj_point_light, inter.normal);     //BRDF

        L_dir = light_sample.emit * f_r * dotProduct(obj_point_light, inter.normal) * dotProduct(-obj_point_light, light_sample.normal)
            / dotProduct(light_sample.coords - inter.coords, light_sample.coords - inter.coords) / pdf;
    }

    //indirect
    if(get_random_float() < RussianRoulette){
        Vector3f nextDir = inter.m->sample(obj_point_light, inter.normal);
        float pdf_indirct = inter.m->pdf(ray.direction, nextDir, inter.normal);
        Vector3f f_r = inter.m->eval(ray.direction, nextDir, inter.normal);

        L_indir = castRay(Ray(inter.coords, nextDir), depth + 1)
            * f_r * dotProduct(nextDir, inter.normal) / pdf_indirct / RussianRoulette;
    }

    return L_dir + L_indir;
    
}
