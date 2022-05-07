using System;

namespace RayTracer
{
    /// <summary>
    /// Class to represent an (infinite) plane in a scene.
    /// </summary>
    public class Sphere : SceneEntity
    {
        private Vector3 center;
        private double radius;
        private Material material;

        /// <summary>
        /// Construct a sphere given its center point and a radius.
        /// </summary>
        /// <param name="center">Center of the sphere</param>
        /// <param name="radius">Radius of the spher</param>
        /// <param name="material">Material assigned to the sphere</param>
        public Sphere(Vector3 center, double radius, Material material)
        {
            this.center = center;
            this.radius = radius;
            this.material = material;
        }

        /// <summary>
        /// Determine if a ray intersects with the sphere, and if so, return hit data.
        /// </summary>
        /// <param name="ray">Ray to check</param>
        /// <returns>Hit data (or null if no intersection)</returns>
        public RayHit Intersect(Ray ray)
        {
            // Write your code here...
            Vector3 oc = ray.Origin-this.center;
            double a = ray.Direction.Dot(ray.Direction);
            double b = 2.0*oc.Dot(ray.Direction);
            double c = oc.Dot(oc) - this.radius*this.radius;
            double disc = b*b - 4*a*c;
            if (disc < 0.0) return null;

            double numerator = -b - Math.Sqrt(disc);
            if (numerator > 0.0){
                Vector3 norm = (ray.Origin + numerator/(2.0*a)*ray.Direction - this.center)/this.radius;
                return new RayHit(ray.Origin + numerator/(2.0*a)*ray.Direction, norm.Normalized(), ray.Direction, this.material);
            }
            numerator = -b + Math.Sqrt(disc);
            if (numerator > 0.0){
                Vector3 norm = (ray.Origin + numerator/(2.0*a)*ray.Direction - this.center)/this.radius;
                return new RayHit(ray.Origin + numerator/(2.0*a)*ray.Direction, norm.Normalized(), ray.Direction, this.material);
            }
            return null;
        }

        /// <summary>
        /// The material of the sphere.
        /// </summary>
        public Material Material { get { return this.material; } }
    }

}
