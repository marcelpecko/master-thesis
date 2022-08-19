using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK;
using MathNet.Numerics.Distributions; // aby sme mohli pouzivat gaussa



namespace ProceduralTerrainGeneration
{
    class Terrain
    {
        public float[,] HeightmapAsArray;
        public float[,] OldHeightmapAsArray;
        public List<Vector3> Coordinates;
        public List<int> Indices;
        public int Iteration;
        public double Hurst;
        public Vector3 Scale;
        public Vector3 OldScale = new Vector3(1, 1, 1);

        // Initialization of a terrain from array
        public Terrain(float[,] _heightmapAsArray, int _iteration, double hurst, Vector3 scale)
        {
            HeightmapAsArray = _heightmapAsArray;
            OldHeightmapAsArray = new float [HeightmapAsArray.GetLength(0),HeightmapAsArray.GetLength(1)];
            for (int i = 0; i < HeightmapAsArray.GetLength(0); i++)
            {
                for (int j = 0; j < HeightmapAsArray.GetLength(1); j++)
                {
                    OldHeightmapAsArray[i, j] = HeightmapAsArray[i, j];
                }
            }
            Iteration = _iteration;
            Hurst = hurst;
            Scale = scale;
            Coordinates = GetCoordinatesFromHeightmap();
            Indices = GetIndices();
            
        }

        private List<Vector3> GetCoordinatesFromHeightmap()
        {            
            int width = HeightmapAsArray.GetLength(1);
            int height = HeightmapAsArray.GetLength(0);
            List<Vector3> coordinatesList = new List<Vector3>();
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                {                           
                    coordinatesList.Add(new Vector3((-2.5f + 5.0f * i / (height - 1))*Scale.X, (-(width - 1) * 2.5f / (height - 1) + 5.0f * j / (height - 1))*Scale.Y, (float)HeightmapAsArray[i, j]*Scale.Z));
                }

            return coordinatesList;
        }

        private List<int> GetIndices()
        {
            int width = HeightmapAsArray.GetLength(1);
            int height = HeightmapAsArray.GetLength(0);
            List<int> indexesList = new List<int>();

            for (int i = 0; i < height - 1; i++)
                for (int j = 0; j < width - 1; j++)
                {
                    // triangles with counterclockwise indexes
                    indexesList.Add(i * width + j);    // 0
                    indexesList.Add((i + 1) * width + j + 1);  //5
                    indexesList.Add(i * width + j + 1);   // 1

                    indexesList.Add(i * width + j);   //0
                    indexesList.Add((i + 1) * width + j);   // 4
                    indexesList.Add((i + 1) * width + j + 1);  //5
                }
            return indexesList;
        }

        public void RecomputeCoordinates()
        {
            Vector3 actualScale = new Vector3(Scale.X/OldScale.X, Scale.Y / OldScale.Y, Scale.Z / OldScale.Z);           
            if (OldScale.X == 0 || OldScale.Y == 0 || OldScale.Z == 0)
            {
                Coordinates.Clear();         

                for (int i = 0; i < HeightmapAsArray.GetLength(0); i++)
                {
                    for (int j = 0; j < HeightmapAsArray.GetLength(1); j++) { 
                        HeightmapAsArray[i, j] = OldHeightmapAsArray[i, j];
                    }
                }                                         

                Coordinates = GetCoordinatesFromHeightmap();
                OldScale = Scale;
                return;
            }
                       
            for (int i = 0; i < Coordinates.Count(); i++)
            {
                Coordinates[i] = Coordinates[i] * actualScale;
            }
            
            for (int i = 0; i < HeightmapAsArray.GetLength(0); i++)
            {
                for (int j = 0; j < HeightmapAsArray.GetLength(1); j++)
                {
                    HeightmapAsArray[i, j] = OldHeightmapAsArray[i, j] * Scale.Z;                   
                }
            }            
            OldScale = Scale;
        }               

        public Terrain DiamondSquareAlgorithm()
        {
            return new Terrain(SquareStep(DiamondStep()), Iteration, Hurst, Scale);
        }

        private float[,] DiamondStep()
        {
            double mean = 0;
            double stdDev = Math.Sqrt(Math.Pow(0.5, Iteration * 2 * Hurst)); // tu asi este ide odmocnina
            MathNet.Numerics.Distributions.Normal normalDist = new Normal(mean, stdDev);

            int oldWidth = HeightmapAsArray.GetLength(1);
            int oldHeight = HeightmapAsArray.GetLength(0);
            if (oldWidth < 2 || oldHeight < 2) return HeightmapAsArray;

            int newWidth = 2 * HeightmapAsArray.GetLength(1) - 1;
            int newHeight = 2 * HeightmapAsArray.GetLength(0) - 1;
            float[,] newHeightmap = new float[newHeight, newWidth];

            for (int i = 0; i < oldHeight; i++)
                for (int j = 0; j < oldWidth; j++)
                {
                    newHeightmap[2 * i, 2 * j] = HeightmapAsArray[i, j];
                    if (j != (oldWidth - 1) && i != (oldHeight - 1))
                    {                        
                        float randomGaussianValue = (float)normalDist.Sample();
                        newHeightmap[2 * i + 1, 2 * j + 1] = (float)(HeightmapAsArray[i, j] + HeightmapAsArray[i + 1, j] + HeightmapAsArray[i, j + 1] + HeightmapAsArray[i + 1, j + 1]) / 4 + randomGaussianValue;
                    }
                }
            return newHeightmap;
        }

        private float[,] SquareStep(float[,] heightmapAsArray)
        {
            double mean = 0;
            double stdDev = Math.Sqrt(Math.Pow(0.5, Iteration * 2 * Hurst)); // tu asi este ide odmocnina
            MathNet.Numerics.Distributions.Normal normalDist = new Normal(mean, stdDev);

            int width = heightmapAsArray.GetLength(1);
            int height = heightmapAsArray.GetLength(0);
            int j = 1;
            for (int i = 0; i < height; i++)
            {
                for (; j < width; j = j + 2)
                {
                    float randomGaussianValue = (float)normalDist.Sample();

                    if (i == 0)
                    {
                        heightmapAsArray[i, j] = (float)(heightmapAsArray[i + 1, j] + heightmapAsArray[i, j + 1] + heightmapAsArray[i, j - 1]) / 3 + randomGaussianValue; 
                    }
                    else if (i == height - 1)
                    {                       
                        heightmapAsArray[i, j] = (float)(heightmapAsArray[i - 1, j] + heightmapAsArray[i, j + 1] + heightmapAsArray[i, j - 1]) / 3 + randomGaussianValue;
                    }
                    else if (j == 0)
                    {
                        heightmapAsArray[i, j] = (float)(heightmapAsArray[i - 1, j] + heightmapAsArray[i + 1, j] + heightmapAsArray[i, j + 1]) / 3 + randomGaussianValue; 
                    }
                    else if (j == width - 1)
                    {
                        heightmapAsArray[i, j] = (float)(heightmapAsArray[i - 1, j] + heightmapAsArray[i + 1, j] + heightmapAsArray[i, j - 1]) / 3 + randomGaussianValue; 
                    }
                    else
                    {
                        heightmapAsArray[i, j] = (float)(heightmapAsArray[i - 1, j] + heightmapAsArray[i + 1, j] + heightmapAsArray[i, j + 1] + heightmapAsArray[i, j - 1]) / 4 + randomGaussianValue;
                    }
                }
                j = (i % 2) != 0 ? 1 : 0;
            }
            return heightmapAsArray;
        }

        public float[,] ReduceSampling()
        {
            int oldWidth = HeightmapAsArray.GetLength(1);
            int oldHeight = HeightmapAsArray.GetLength(0);

            int newWidth = (int)Math.Ceiling((float)oldWidth / 2);
            int newHeight = (int)Math.Ceiling((float)oldHeight / 2);
            if (newWidth < 2) newWidth = oldWidth;
            if (newHeight < 2) newHeight = oldHeight;
            float[,] newHeightmap = new float[newHeight, newWidth];

            for (int i = 0; i < oldHeight; i = i + 2)
                for (int j = 0; j < oldWidth; j = j + 2)
                {
                    newHeightmap[i / 2, j / 2] = HeightmapAsArray[i, j];
                }
            return newHeightmap;
        }
    }
}
