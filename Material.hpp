//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"

enum MaterialType { DIFFUSE, Microfacet};

class Material{
private:

    // Compute reflection direction
    Vector3f reflect(const Vector3f &I, const Vector3f &N) const
    {
        return I - 2 * dotProduct(I, N) * N;
    }

    // Compute refraction direction using Snell's law
    //
    // We need to handle with care the two possible situations:
    //
    //    - When the ray is inside the object
    //
    //    - When the ray is outside.
    //
    // If the ray is outside, you need to make cosi positive cosi = -N.I
    //
    // If the ray is inside, you need to invert the refractive indices and negate the normal N
    Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        Vector3f n = N;
        if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }

    // Compute Fresnel equation
    //
    // \param I is the incident view direction
    //
    // \param N is the normal at the intersection point
    //
    // \param ior is the material refractive index
    //
    // \param[out] kr is the amount of light reflected

    /*  //完整的菲涅尔方程，fresnel equation
    void fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        if (cosi > 0) {  std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            kr = 1;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            kr = (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is given by:
        // kt = 1 - kr;
    }
    */

    //简化版菲涅尔方程， simplified fresnel equation
    float fresnelSchlick(const float &n1, Vector3f L, Vector3f H){
        float R0 = pow((n1-1)/(n1+1), 2);
        return R0 + (1.0 - R0) * pow(1.0 - std::max(dotProduct(L, H), 0.0f), 5);
    }


    Vector3f toWorld(const Vector3f &a, const Vector3f &N){
        Vector3f B, C;
        if (std::fabs(N.x) > std::fabs(N.y)){
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            C = Vector3f(N.z * invLen, 0.0f, -N.x *invLen);
        }
        else {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C = Vector3f(0.0f, N.z * invLen, -N.y *invLen);
        }
        B = crossProduct(C, N);
        return a.x * B + a.y * C + a.z * N;
    }

private:
    float DistributionGGX(Vector3f N, Vector3f H, float roughness)
    {
        float a = roughness * roughness;
        float a2 = a * a;
        float NdotH = std::max(dotProduct(N, H), 0.0f);
        float NdotH2 = NdotH * NdotH;

        float nom = a2;
        float denom = (NdotH2 * (a2 - 1.0) + 1.0);
        denom = M_PI * denom * denom;

        return nom / std::max(denom, 0.0000001f); // prevent divide by zero for roughness=0.0 and NdotH=1.0
    }

    float GeometrySchlickGGX(float NdotV, float k)
    {
        float nom = NdotV;
        float denom = NdotV * (1.0 - k) + k;
        return nom / denom;
    }

    float GeometrySmith(Vector3f N, Vector3f V, Vector3f L, float roughness)
    {
        float r = (roughness + 1.0);
        float k = (r * r) / 8.0;
        float NdotV = std::max(dotProduct(N, V), 0.0f);
        float NdotL = std::max(dotProduct(N, L), 0.0f);
        float ggx2 = GeometrySchlickGGX(NdotV, k);
        float ggx1 = GeometrySchlickGGX(NdotL, k);

        return ggx1 * ggx2;
    }


public:
    MaterialType m_type;
    //Vector3f m_color;
    Vector3f m_emission;
    float ior;
    Vector3f Kd, Ks;
    float specularExponent;
    //Texture tex;

    inline Material(MaterialType t=DIFFUSE, Vector3f e=Vector3f(0,0,0));
    inline MaterialType getType();
    //inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();
    inline bool hasEmission();

    // sample a ray by Material properties
    inline Vector3f sample(const Vector3f &wi, const Vector3f &N);
    // given a ray, calculate the PdF of this ray
    inline float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
    // given a ray, calculate the contribution of this ray
    inline Vector3f eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);

};

Material::Material(MaterialType t, Vector3f e){
    m_type = t;
    //m_color = c;
    m_emission = e;
}

MaterialType Material::getType(){return m_type;}
///Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() {return m_emission;}
bool Material::hasEmission() {
    if (m_emission.norm() > EPSILON) return true;
    else return false;
}

Vector3f Material::getColorAt(double u, double v) {
    return Vector3f();
}


Vector3f Material::sample(const Vector3f &wi, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r*std::cos(phi), r*std::sin(phi), z);
            return toWorld(localRay, N);
            
            break;
        }
        case Microfacet:
        {
            // uniform sample on the hemisphere在半球上均匀采样
            float x_1 = get_random_float(), x_2 = get_random_float();
            //z∈[0,1]，是随机半球方向的z轴向量
            float z = std::fabs(1.0f - 2.0f * x_1);
            //r是半球半径随机向量以法线为旋转轴的半径
            //phi是r沿法线旋转轴的旋转角度
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;//phi∈[0,2*pi]
            Vector3f localRay(r*std::cos(phi), r*std::sin(phi), z);//半球面上随机的光线的弹射方向
            return toWorld(localRay, N);//转换到世界坐标
            
            break;
        }
    }
}

float Material::pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample probability 1 / (2 * PI)
            if (dotProduct(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
        case Microfacet:
        {
            // uniform sample probability 1 / (2 * PI)
            if (dotProduct(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
    }
}

Vector3f Material::eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // calculate the contribution of diffuse   model
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f) {
                Vector3f diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return Vector3f(0.0f);
            break;
        }
        case Microfacet://微表面材质的BRDF
        {
            // Disney PBR 方案
            float cosalpha = dotProduct(N, wo);
            if (cosalpha > 0.0f) {
                float roughness = 0.40;

                Vector3f V = -wi;
                Vector3f L = wo;
                Vector3f H = normalize(V + L);

                // 计算 distribution of normals: D
                float D = DistributionGGX(N, H, roughness);

                // 计算 shadowing masking term: G
                float G = GeometrySmith(N, V, L, roughness);

                // 计算 fresnel 系数: F
                float F;
                float etat = 1.85;
                //fresnel(wi, N, etat, F);
                F = fresnelSchlick(etat, L, H);

                Vector3f nominator = D * G * F;
                float denominator = 4 * std::max(dotProduct(N, V), 0.0f) * std::max(dotProduct(N, L), 0.0f);
                Vector3f specular = nominator / std::max(denominator, 0.001f);

                // 能量守恒
                float ks_ = F;//反射比率
                float kd_ = 1.0f - ks_;//折射比率

                Vector3f diffuse = 1.0f / M_PI;

                // 因为在 specular 项里已经考虑了反射部分的比例：F。所以反射部分不需要再乘以 ks_ 
                //Ks为镜面反射项，Kd为漫反射项。
                return Ks * specular + kd_ * Kd * diffuse;
            }
            else
                return Vector3f(0.0f);
            break;
        }
    }
}


#endif //RAYTRACING_MATERIAL_H
