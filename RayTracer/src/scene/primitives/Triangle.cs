using System;

namespace RayTracer
{
    /// <summary>
    /// Class to represent a triangle in a scene represented by three vertices.
    /// </summary>
    public class Triangle : SceneEntity
    {
        private Vector3 v0, v1, v2;
        private Material material;

        /// <summary>
        /// Construct a triangle object given three vertices.
        /// </summary>
        /// <param name="v0">First vertex position</param>
        /// <param name="v1">Second vertex position</param>
        /// <param name="v2">Third vertex position</param>
        /// <param name="material">Material assigned to the triangle</param>
        public Triangle(Vector3 v0, Vector3 v1, Vector3 v2, Material material)
        {
            this.v0 = v0;
            this.v1 = v1;
            this.v2 = v2;
            this.material = material;
        }

        /// <summary>
        /// Determine if a ray intersects with the triangle, and if so, return hit data.
        /// </summary>
        /// <param name="ray">Ray to check</param>
        /// <returns>Hit data (or null if no intersection)</returns>
        public RayHit Intersect(Ray ray)
        {
            Vector3 norm = (v1-v0).Cross(v2-v0);
            if (Math.Abs(ray.Direction.Dot(norm)) >= 0.001){
                double d = norm.Dot(v0-ray.Origin);
                double t = d / norm.Dot(ray.Direction);
                Vector3 p = ray.Origin + t*ray.Direction;

                Vector3 a = (v1-v0).Cross(p-v0);
                Vector3 b = (v2-v1).Cross(p-v1);
                Vector3 c = (v0-v2).Cross(p-v2);

                if ((norm.Dot(a) < 0) || (norm.Dot(b) < 0) || (norm.Dot(c) < 0)){
                    return null;
                }
                else{
                    if(t >= 0.01){
                        RayHit ob = new RayHit(p, norm.Normalized(), ray.Direction, this.material);
                        return ob;
                    }
                }
            }
            return null;
        }

        /// <summary>
        /// The material of the triangle.
        /// </summary>
        public Material Material { get { return this.material; } }
    }

}
