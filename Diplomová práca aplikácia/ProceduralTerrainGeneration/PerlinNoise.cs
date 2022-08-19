using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using OpenTK;


namespace ProceduralTerrainGeneration
{
    class PerlinNoise
    {
        public float[,] HeightmapAsArray;

        public PerlinNoise(int width, int height, float scale, int seed) {
            HeightmapAsArray = new float[height, width];

            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    //HeightmapAsArray[i, j] = (double)perlin((float)i /(height-1), (float)j/(width-1));
                    HeightmapAsArray[i, j] = perlin(i * scale, j * scale, seed);
        }


        float interpolate(float a0, float a1, float w)
        {
            /* // You may want clamping by inserting:
             * if (0.0 > w) return a0;
             * if (1.0 < w) return a1;
             */
            return (a1 - a0) * w + a0;
            /* // Use this cubic interpolation [[Smoothstep]] instead, for a smooth appearance:
             * return (a1 - a0) * (3.0 - w * 2.0) * w * w + a0;
             *
             * // Use [[Smootherstep]] for an even smoother result with a second derivative equal to zero on boundaries:
             * return (a1 - a0) * ((w * (w * 6.0 - 15.0) + 10.0) * w * w * w) + a0;
             */
        }

        Vector2 randomGradient(int ix, int iy, float randomFloat)
        {
            // No precomputed gradients mean this works for any number of grid coordinates
            const int w = 8 * sizeof(int);
            const int s = w / 2; // rotation width
            long a = ix, b = iy;
            a *= 3284157443; b ^= a << s | a >> w - s;
            b *= 1911520717; a ^= b << s | b >> w - s;
            a *= 2048419325;
                                  
            float random = (float)(a* randomFloat * (3.14159265 / ~(~0u >> 1))); // in [0, 2*Pi]
            Vector2 v;
            
            v.X = (float)Math.Cos(random); 
            v.Y = (float)Math.Sin(random);
            return v;
        }

        float dotGridGradient(int ix, int iy, float x, float y, float randomFloat)
        {
            // Get gradient from integer coordinates
            Vector2 gradient = randomGradient(ix, iy, randomFloat);

            // Compute the distance vector
            float dx = x - (float)ix;
            float dy = y - (float)iy;

            // Compute the dot-product
            return (dx * gradient.X + dy * gradient.Y);
        }

        float perlin(float x, float y, int seed)
        {
            var random = new Random(seed);
            float randomFloat = (float)random.NextDouble(); 

            // Determine grid cell coordinates
            int x0 = (int)Math.Floor(x);
            int x1 = x0 + 1;
            int y0 = (int)Math.Floor(y);
            int y1 = y0 + 1;

            // Determine interpolation weights
            // Could also use higher order polynomial/s-curve here
            float sx = x - (float)x0;
            float sy = y - (float)y0;

            // Interpolate between grid point gradients
            float n0, n1, ix0, ix1, value;

            n0 = dotGridGradient(x0, y0, x, y, randomFloat);
            n1 = dotGridGradient(x1, y0, x, y, randomFloat);
            ix0 = interpolate(n0, n1, sx);

            n0 = dotGridGradient(x0, y1, x, y, randomFloat);
            n1 = dotGridGradient(x1, y1, x, y, randomFloat);
            ix1 = interpolate(n0, n1, sx);

            value = interpolate(ix0, ix1, sy);
            return value;
        }      

    }
}
