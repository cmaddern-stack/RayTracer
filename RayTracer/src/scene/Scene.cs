using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace RayTracer
{
    /// <summary>
    /// Class to represent a ray traced scene, including the objects,
    /// light sources, and associated rendering logic.
    /// </summary>
    public class Scene
    {
        private SceneOptions options;
        private ISet<SceneEntity> entities;
        private ISet<PointLight> lights;
        private Color noLight = new Color(0,0,0);

        /// <summary>
        /// Construct a new scene with provided options.
        /// </summary>
        /// <param name="options">Options data</param>
        public Scene(SceneOptions options = new SceneOptions())
        {
            this.options = options;
            this.entities = new HashSet<SceneEntity>();
            this.lights = new HashSet<PointLight>();
        }

        /// <summary>
        /// Add an entity to the scene that should be rendered.
        /// </summary>
        /// <param name="entity">Entity object</param>
        public void AddEntity(SceneEntity entity)
        {
            this.entities.Add(entity);
        }

        /// <summary>
        /// Add a point light to the scene that should be computed.
        /// </summary>
        /// <param name="light">Light structure</param>
        public void AddPointLight(PointLight light)
        {
            this.lights.Add(light);
        }

        /// <summary>
        /// Determine the Color of a diffusion object. Use the value of the rayhit 
        /// and the entity colided with to determine color.
        /// </summary>
        /// <param name="hit">Point of ray intersection</param>
        /// <param name="entity">Entity collided with</param>
        private Color HandleDiffusion(RayHit hit, SceneEntity entity)
        {
            // If no light rays intersect then return black
            if(this.lights.Count == 0){
                return noLight;
            }
            else{
                Color c_l = noLight;
                int first = 0;
                foreach (PointLight light in this.lights){

                    // Determine if the path between a light and the object is blocked
                    double distance = Math.Abs((light.Position-hit.Position).Length());
                    bool isBlocked = false;
                    Ray check = new Ray(hit.Position, (light.Position-hit.Position).Normalized());
                    foreach(SceneEntity block in this.entities){
                        RayHit block_h = block.Intersect(check);
                        if(block_h != null){
                            double d2 = Math.Abs((block_h.Position - hit.Position).Length());
                            if((d2 < distance && d2 > 0.0001)){
                                isBlocked = true;
                                break;
                            }
                        }
                    }
                    // If light is blocked then dont contribute colour of light
                    if (isBlocked) continue;

                    // Add contribution of given light to diffusion output
                    if (first == 0){
                        double val = hit.Normal.Dot((light.Position - hit.Position).Normalized());
                        if (val < 0) val = 0;
                        c_l = entity.Material.Color*light.Color*val;
                        first = 1;
                    } 
                    else{
                        double val = hit.Normal.Dot((light.Position - hit.Position).Normalized());
                        if (val < 0) val = 0;
                        c_l += entity.Material.Color*light.Color*val;
                    }
                }
                return c_l;
            }
        }


        /// <summary>
        /// Determine the Color of a reflection object. Use the value of the rayhit 
        /// and the entity colided with to determine color.
        /// </summary>
        /// <param name="hit">Point of ray intersection</param>
        /// <param name="ray">Original ray object</param>
        /// <param name="depth">Number of iterations to explore</param>
        /// <param name="SceneEntity">Entity collided with</param>
        private Color HandleReflection(RayHit hit, Ray ray, int depth, SceneEntity entity)
        {
            Vector3 reflect_dir = ray.Direction - 2*(ray.Direction.Dot(hit.Normal)*hit.Normal);
            Ray reflection_ray = new Ray(hit.Position+reflect_dir*0.00001, reflect_dir);

            return HandleCollisions(reflection_ray, depth);
        }


        /// <summary>
        /// Determine the Color of a refraction object. Use the value of the rayhit 
        /// the original ray and the entities refractive index to determine color.
        /// </summary>
        /// <param name="hit">Point of ray intersection</param>
        /// <param name="ray">Original ray object</param>
        /// <param name="depth">Number of iterations to explore</param>
        /// <param name="ior">Refractive index of entity</param>
        private Color HandleRefraction(RayHit hit, Ray ray, int depth, double ior)
        {
            double cosi = Math.Clamp(ray.Direction.Dot(hit.Normal), -1, 1);
            
            double etai = 1, etat = ior;
            Vector3 n = hit.Normal;

            // Check if index is outside of the surface
            if (cosi < 0){
                cosi = -cosi;
            }
            else{
                // Swap refraction indicies
                etai = etat;
                etat = 1; 
                n = -hit.Normal;
            }
            double eta = etai/etat;
            double k = 1 - eta*eta*(1-cosi*cosi);
            Vector3 refract_dir;
            if (k < 0){
                return new Color(0,0,0);
            }
            else{
                refract_dir = eta*ray.Direction+(eta*cosi-Math.Sqrt(k))*n;
            }

            Ray refraction_ray = new Ray(hit.Position+refract_dir*0.00001, refract_dir);
            
            
            return HandleCollisions(refraction_ray, depth);
        }


        /// <summary>
        /// Determine the type of oject a ray collides with and handle possible
        /// reflection / refraction etc.
        /// </summary>
        /// <param name="resultant_ray">Ray generated from previous collision</param>
        /// <param name="depth">Number of iterations to explore</param>
        public Color HandleCollisions(Ray resultant_ray, int depth){
            double distance = double.PositiveInfinity;
            SceneEntity closest_entity = null;
            RayHit closest_hit = null;
            Color ref_col = noLight;

            
            // Determine the closest entity of which the ray intersects with
            foreach(SceneEntity entity_ref in this.entities){
                RayHit hit_ray = entity_ref.Intersect(resultant_ray);
                if (hit_ray != null){
                    double ref_dis = (resultant_ray.Origin - hit_ray.Position).Length();
                    if( ref_dis < distance){
                        closest_entity = entity_ref;
                        closest_hit = hit_ray;
                        distance = ref_dis;
                    }
                }
            }

            if(closest_entity != null){
                // Handle collision with a diffuse object
                if (closest_entity.Material.Type is Material.MaterialType.Diffuse){
                    // Check if implement ambient lighting
                    if(options.AmbientLightingEnabled)
                        ref_col = castRay(resultant_ray.Origin, resultant_ray.Direction, 3, closest_entity, closest_hit);
                    else
                        ref_col = HandleDiffusion(closest_hit, closest_entity);
                    
                }
                // Handle collision with a reflective object
                else if( closest_entity.Material.Type is Material.MaterialType.Reflective){
                    if (depth >= 0){
                        ref_col = HandleReflection(closest_hit, resultant_ray, depth-1, closest_entity);
                    }
                }
                // Handle collision with a refractive object
                else if( closest_entity.Material.Type is Material.MaterialType.Refractive){
                    if (depth >= 0){
                        ref_col = HandleFresnel(resultant_ray, closest_hit, closest_entity, depth-1); 
                    }
                }
                else{
                    ref_col= closest_entity.Material.Color;
                }
            }

            return ref_col;
        }

        /// <summary>
        /// Determine the ratio of reflection to refraction on a medium using the fresnel equations
        /// </summary>
        /// <param name="ray">Original ray</param>
        /// <param name="hit">Point of ray intersection</param>
        /// <param name="entity">Entity collided with</param>
        /// <param name="depth">Number of iterations to explore</param>
        public Color HandleFresnel(Ray ray, RayHit hit, SceneEntity entity, int depth){

            // Determine refraction and reflection contribution
            Color refract = HandleRefraction(hit, ray, depth, entity.Material.RefractiveIndex);
            Color reflect = HandleReflection(hit, ray, depth, entity);
            
            double kr, kt;

            double cosi = Math.Clamp(ray.Direction.Dot(hit.Normal), -1, 1);
            double etai = 1, etat = entity.Material.RefractiveIndex;
            if (cosi > 0){
                etai = etat;
                etat = 1;
            }
            // Compute Snell's law for sint
            double sint = (etai / etat) * Math.Sqrt(Math.Max(0.0, 1-cosi*cosi));
            //Check for total internal reflection
            if(sint >= 1) {
                kr = 1;
            }
            else{
                double cost = Math.Sqrt(Math.Max(0.0, 1 - sint*sint));
                cosi = Math.Abs(cosi);
                double Rs = ((etat*cosi)-(etai*cost)) / ((etat*cosi)+(etai*cost));
                double Rp = ((etai*cosi)-(etat*cost)) / ((etai*cosi)+(etat*cost));
                kr = (Rs*Rs + Rp * Rp) / 2;
            }
            // Determine ratio of reflection to refraction
            kt = 1 - kr;

            return refract*kt+reflect*kr;
        }

        public Color castRay(Vector3 orig, Vector3 dir, int depth, SceneEntity entity, RayHit hit){
            if (depth <= 0){
                return noLight; 
            } 

            // Direct lighting contribution
            Color directLightContrib = HandleDiffusion(hit, entity);

            double indirect_r = 0.0;
            double indirect_g = 0.0;
            double indirect_b = 0.0;


            Vector3 Nt, Nb;

            // Create coordinate system along shading normal
            Tuple<Vector3,Vector3> co = createCoordinateSystem(hit.Normal);

            Nt = co.Item1;
            Nb = co.Item2;

            double pdf = 1/(2*Math.PI);

            int num = 16;
            for (int n = 0; n < num; n++){

                // Sample rays in the hemisphere
                Random rnd = new Random();
                double r1 = rnd.NextDouble();
                double r2 = rnd.NextDouble();
                Vector3 sample = uniformSampleHemisphere(r1, r2);

                Vector3 sampleWorld = new Vector3(
                    sample.X*Nb.X + sample.Y*hit.Normal.X + sample.Z * Nt.X,
                    sample.X*Nb.Y + sample.Y*hit.Normal.Y + sample.Z * Nt.Y,
                    sample.X*Nb.Z + sample.Y*hit.Normal.Z + sample.Z * Nt.Z
                );

                // Determine contributing entity
                Ray sampleRay = new Ray(hit.Position + sampleWorld*0.00001, sampleWorld);
                double distance = double.PositiveInfinity;
                SceneEntity closest = null;
                RayHit closest_hit = null;
                foreach (SceneEntity entity_hit in this.entities){
                    RayHit sample_hit = entity_hit.Intersect(sampleRay);
                    if (sample_hit != null){
                        double ref_dis = (sampleRay.Origin - sample_hit.Position).Length();
                        if (ref_dis < distance){
                            closest = entity_hit;
                            closest_hit = sample_hit;
                            distance = ref_dis;
                        }
                    }
                }
                // Calculate indirect lighting contribution from sample
                if (closest != null && closest.Material.Type is Material.MaterialType.Diffuse){
                    Color found = HandleDiffusion(closest_hit, closest);
                    indirect_r += found.R*r1/pdf ;
                    indirect_g += found.G*r1/pdf;
                    indirect_b += found.B*r1/pdf;
                }
            }
            
            indirect_r /= num;
            indirect_g /= num;
            indirect_b /= num;

            indirect_r += directLightContrib.R;
            indirect_g += directLightContrib.G;
            indirect_b += directLightContrib.B;

            indirect_r /= Math.PI;
            indirect_g /= Math.PI;
            indirect_b /= Math.PI;

            Color indirectLightContrib = new Color(indirect_r, indirect_g, indirect_b);

            return indirectLightContrib * entity.Material.Color;
        }


        /// <summary>
        /// Create a custom coordinate system
        /// </summary>
        /// <param name="N">entity shader normal</param>
        public Tuple<Vector3, Vector3> createCoordinateSystem(Vector3 N){
            Vector3 Nt, Nb;
            if(Math.Abs(N.X) > Math.Abs(N.Y)){
                Nt = new Vector3(N.Z, 0, -N.X) / Math.Sqrt(N.X*N.X+N.Z*N.Z);
            }
            else{
                Nt = new Vector3(0, -N.Z, N.Y) / Math.Sqrt(N.Y*N.Y+N.Z*N.Z);
            }
            Nb = N.Cross(Nt);
            return Tuple.Create(Nt, Nb);
        }

        /// <summary>
        /// Compute uniformly distrubted direction over a hemisphere
        /// </summary>
        /// <param name="r1">Random variable betwwen 0-1</param>
        /// <param name="r2">Random variable betwwen 0-1</param>
        public Vector3 uniformSampleHemisphere(double r1, double r2){
            double sinTheta = Math.Sqrt(1-r1*r1);
            double phi = 2*Math.PI*r2;
            double x = sinTheta * Math.Cos(phi);
            double z = sinTheta * Math.Sin(phi);
            return new Vector3(x, r1, z);
        }



        


        /// <summary>
        /// Render the scene to an output image. This is where the bulk
        /// of your ray tracing logic should go... though you may wish to
        /// break it down into multiple functions as it gets more complex!
        /// </summary>
        /// <param name="outputImage">Image to store render output</param>
        public void Render(Image outputImage)
        {
            // Begin writing your code here...
            var watch = System.Diagnostics.Stopwatch.StartNew();
            int res_val = options.AAMultiplier;

            // Use parallel for to run on all cores simultaneously to boost speed of render
            Parallel.For(0, outputImage.Height, y =>
            {
                for(int x = 0; x < outputImage.Width; x++){
                    double aa_r = 0, aa_g = 0, aa_b = 0;
                    for(int x_res = 0; x_res < res_val; x_res++){
                        for(int y_res = 0; y_res < res_val; y_res++){


                            // Normalize the pixel space
                            double pixel_loc_x = (res_val*x+x_res + 0.5 ) / ((double)outputImage.Width*res_val);
                            double pixel_loc_y = (res_val*y+y_res + 0.5) / ((double)outputImage.Height*res_val);
                            double pixel_loc_z = 1.0f;

                            double x_pos = (pixel_loc_x * 2) - 1;
                            double y_pos = 1 - (pixel_loc_y * 2);

                            // Set FOV
                            x_pos = x_pos * Math.Tan(Math.PI / 6);
                            y_pos = y_pos * (Math.Tan(Math.PI / 6) / ((double)outputImage.Width / (double)outputImage.Height));

                            Ray ray = new Ray(new Vector3(0, 0, 0), (new Vector3(x_pos, y_pos, pixel_loc_z).Normalized()));
                            
                            double z_in = double.PositiveInfinity;
                            Color final_color = noLight;


                            // Find first entity ray collides with
                            SceneEntity closest_entity = null;
                            RayHit closest_hit = null;

                            foreach (SceneEntity entity in this.entities){
                                RayHit hit = entity.Intersect(ray);
                                if (hit != null){
                                    if(z_in >= hit.Position.Z){
                                        closest_entity = entity;
                                        closest_hit = hit;
                                        z_in = hit.Position.Z;
                                    }
                                }
                            }
                            if(closest_entity != null){
                                // Handle color depending on collision material type
                                if (closest_entity.Material.Type is Material.MaterialType.Diffuse){
                                    if(options.AmbientLightingEnabled)
                                        final_color = castRay(ray.Origin, ray.Direction, 3, closest_entity, closest_hit);
                                    else
                                        final_color = HandleDiffusion(closest_hit, closest_entity);
                                    outputImage.SetPixel(x,y, final_color);
                                }
                                else if (closest_entity.Material.Type is Material.MaterialType.Reflective){
                                    final_color = HandleReflection(closest_hit, ray, 3, closest_entity);
                                    outputImage.SetPixel(x,y, final_color);
                                }
                                else if(closest_entity.Material.Type is Material.MaterialType.Refractive){
                                    final_color = HandleFresnel(ray, closest_hit, closest_entity, 3);
                                    outputImage.SetPixel(x, y, final_color);
                                }
                                else if(closest_entity.Material.Type is Material.MaterialType.Emissive){
                                    final_color = closest_entity.Material.Color;
                                    outputImage.SetPixel(x, y, final_color);
                                }
                                else{
                                    final_color = closest_entity.Material.Color;
                                    outputImage.SetPixel(x,y,closest_entity.Material.Color);
                                }
                            }
                            
                            // Add color for averaging when using greater aa values
                            aa_r += final_color.R;
                            aa_g += final_color.G;
                            aa_b += final_color.B;
                        }
                    }
                    Color avgColor = new Color(aa_r/(res_val*res_val),aa_g/(res_val*res_val),aa_b/(res_val*res_val));
                    outputImage.SetPixel(x,y,avgColor);
                }
            });
        watch.Stop();
        Console.WriteLine(watch.ElapsedMilliseconds);
        }
    }
}
