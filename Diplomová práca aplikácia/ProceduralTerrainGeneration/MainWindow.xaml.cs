using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Forms;
using System.Windows.Forms.Integration;
using System.Windows.Media.Media3D;
using Microsoft.Win32;
using System.IO;
using System.Drawing;

// OpenTK libraries
using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;
// library for Gauss distribution
using MathNet.Numerics.Distributions;

namespace ProceduralTerrainGeneration
{
    public partial class MainWindow : Window
    {
        // global parameters
        GLControl glControl;
        Terrain MainTerrain = null;

        // cameera settings
        double Dist = new double(), Phi = new double(), Theta = new double(), oPhi = new double(), oTheta = new double(), prevPhi = new double(), prevTheta = new double(), prevDist = new double();

        // mouse settings
        double RightX, RightY;
        bool IsLeftDown, IsRightDown;
        int SelectedPoint;

        // keyboard settings
        bool IsZ = true, IsY = false, IsX = false;
        bool DrawCoordinateAxes = false;
        bool DrawMesh = true;

        // terrain properties
        int NumberOfPointsWidth;
        int NumberOfPointsHeight;
        int NumberOfIterations;

        // water level properties
        bool AddWaterLevel = false;
        int IndexOfMovedPointOfWaterLevel = -1;        
        List<System.Windows.Point> waterLevelPoints;

        // scale properties
        Vector3 Scale;
        float oldZscale = 1;
        float oldWaterLevelValue = 0;

        // colors of elements in application
        readonly float[] WaterColor = { 0.1f, 0.8f, 1.0f };     //RGB (154,154,102)
        readonly float[] TerrainColor = { 0.6f, 0.6f, 0.4f };   //RGB (154,154,102)
        readonly float[] NetColor = { 0.3f, 0.3f, 0.2f, 1.0f }; //RGB (77,77,51)

        // others parameter
        bool PerlinNoiseActive = false;
        bool SimplexNoiseActive = false;
        bool MathFunctionActive = false;
        bool WasChangeOfIteration = false;
        bool sizeSliderDrag = false;

        public MainWindow(string terrainType)
        {
            InitializeComponent();
            float[,] heightmapAsArray = new float[0, 0];
            bool isTerrain = true;
            switch (terrainType)
            {
                case "perlin":
                    isTerrain = false;
                    PerlinNoise perlinNoise = new PerlinNoise((int)sliderWidthHeightmap.Value*4, (int)sliderHeightHeightmap.Value * 4, (float)noiseScaleSlider.Value, (int)noiseSeedSlider.Value); // if we set small width/height values, we would have few noise samples and we would not see the noise
                    heightmapAsArray = ReduceHeightmap(ReduceHeightmap(perlinNoise.HeightmapAsArray));  // we want to work with fewer samples so we can edit the terrain so, we reduce terrain twice
                    terrainType = "perlinov šum";
                    PerlinNoiseActive = true;
                    break;
                case "simplex":
                    isTerrain = false;
                    SimplexNoise simplexNoise = new SimplexNoise((int)sliderWidthHeightmap.Value * 4, (int)sliderHeightHeightmap.Value * 4, (float)noiseScaleSlider.Value, (int)noiseSeedSlider.Value); // if we set small width/height values, we would have few noise samples and we would not see the noise 
                    heightmapAsArray = ReduceHeightmap(ReduceHeightmap(simplexNoise.HeightmapAsArray)); // we want to work with fewer samples so we can edit the terrain so, we reduce terrain twice
                    terrainType = "simplexný šum";
                    SimplexNoiseActive = true;
                    break;
                case "open":
                    isTerrain = false;
                    OpenFile(); return;
                case "math":
                    isTerrain = false;
                    int functionIndex = new Random().Next(functionInCombobox.Items.Count);
                    functionInCombobox.SelectedItem = functionInCombobox.Items[functionIndex];
                    heightmapAsArray = GetHeightmapFromFunction(functionIndex);
                    terrainType = functionInCombobox.Text;
                    MathFunctionActive = true;
                    break;
            }

            if (isTerrain)
            {
                string applicationPath = Environment.CurrentDirectory;
                string pathToHeightmap = applicationPath.Remove(applicationPath.Length - 9) + "Heightmaps/" + terrainType + ".png"; //leave actual path and open file Heightmaps/                
                Bitmap heightmapAsImage = (Bitmap)System.Drawing.Image.FromFile(pathToHeightmap, true);
                heightmapAsArray = GetHeightmapFromImage(heightmapAsImage);
            }
            InitializeTerrainParameters(heightmapAsArray, terrainType);
        }

        private void InitializeTerrainParameters(float[,] mapa, string terrainType)
        {
            IsLeftDown = false;
            IsRightDown = false;

            // setting default scale parameters
            XScaleSlider.Value = 1;
            YScaleSlider.Value = 1;
            ZScaleSlider.Value = 1;
            Scale = new Vector3(1, 1, 1);
            
            // setting default water level parameters
            waterLevelPoints = new List<System.Windows.Point>();
            waterLevelPointsPicker.Children.Clear();
            waterLevelPointsPicker.IsEnabled = false;
            waterLevelSlider.Value = 0;
            waterLevelSlider.IsEnabled = false;
            waterLevelButton.Content = "Pridaj vodnú hladinu";
            waterLevelCutButton.IsEnabled = false;

            // setting default noise parameters
            noiseSeedSlider.IsEnabled = false;
            noiseScaleSlider.IsEnabled = false;

            // setting default functions parameters
            rangeSliderMinX.IsEnabled = false;
            rangeSliderMaxX.IsEnabled = false;
            rangeSliderMinY.IsEnabled = false;
            rangeSliderMaxY.IsEnabled = false;

            if (PerlinNoiseActive || SimplexNoiseActive) {
                noiseSeedSlider.IsEnabled = true;
                noiseScaleSlider.IsEnabled = true;
            }

            if (MathFunctionActive)
            {
                rangeSliderMinX.IsEnabled = true;
                rangeSliderMaxX.IsEnabled = true;
                rangeSliderMinY.IsEnabled = true;
                rangeSliderMaxY.IsEnabled = true;
            }

            WasChangeOfIteration = false;

            NumberOfPointsWidth = mapa.GetLength(1);
            NumberOfPointsHeight = mapa.GetLength(0);
            NumberOfIterations = NumberOfPointsWidth >= NumberOfPointsHeight ? (int)Math.Log(NumberOfPointsWidth, 2) : (int)Math.Log(NumberOfPointsHeight, 2);
                        
            sliderHeightHeightmap.Value = NumberOfPointsHeight;
            sliderWidthHeightmap.Value = NumberOfPointsWidth;

            // add text into statusbar
            terrainTypeNameTextBlock.Text = terrainType.ToUpper();
            terrainInfoTextBlock.Text = "šírka terénu je " + NumberOfPointsWidth + " bodov a dĺžka je " + NumberOfPointsHeight + " bodov";
            terrainInfoTextBlock2.Text = "počet trojuholníkov siete je " + (2 * (NumberOfPointsWidth - 1) * (NumberOfPointsHeight - 1));
            iterationLabel.Content = Convert.ToString(NumberOfIterations);
            MainTerrain = new Terrain(mapa, NumberOfIterations, hurstSlider.Value, Scale);
        }

        //----------------------------------OPENTK AND THINGS NEEDED FOR RENDERING-----------------------------------------------

        private void RenderingAreaSettings(object sender, EventArgs e)
        {
            // OpenTK initialization
            OpenTK.Toolkit.Init();
            var flags = GraphicsContextFlags.Default;
            glControl = new GLControl(new GraphicsMode(32, 24), 2, 0, flags);
            glControl.MakeCurrent();
            glControl.Paint += RenderingSettings;
            glControl.Dock = DockStyle.Fill;
            (sender as WindowsFormsHost).Child = glControl;

            // user controls
            glControl.MouseDown += GLControl_MouseDown;
            glControl.MouseMove += GLControl_MouseMove;
            glControl.MouseUp += GLControl_MouseUp;
            glControl.MouseWheel += GLControl_MouseWheel;

            // shading
            GL.ShadeModel(ShadingModel.Smooth);

            // color of the window
            GL.ClearColor(1.0f, 1.0f, 1.0f, 1.0f);
            GL.ClearDepth(1.0f);

            //enable z-buffering
            GL.Enable(EnableCap.DepthTest);
            GL.DepthFunc(DepthFunction.Lequal);
            GL.Hint(HintTarget.PerspectiveCorrectionHint, HintMode.Nicest);

            GL.Enable(EnableCap.Blend);
            GL.BlendFunc(BlendingFactorSrc.SrcAlpha, BlendingFactorDest.OneMinusSrcAlpha);

            //smoothing
            GL.Enable(EnableCap.LineSmooth);
            GL.Enable(EnableCap.PointSmooth);

            // illumination
            float[] light_ambient = { 0.3f, 0.3f, 0.3f, 1.0f };
            float[] light_diffuse = { 0.4f, 0.4f, 0.4f, 0.0f };
            float[] light_specular = { 0.5f, 0.5f, 0.5f, 1.0f };
            float[] light_position = { 10.0f, 10.0f, 200.0f };
            GL.Light(LightName.Light0, LightParameter.Ambient, light_ambient);
            GL.Light(LightName.Light0, LightParameter.Diffuse, light_diffuse);
            GL.Light(LightName.Light0, LightParameter.Specular, light_specular);
            GL.Light(LightName.Light0, LightParameter.ConstantAttenuation, 1.0f);
            GL.Light(LightName.Light0, LightParameter.Position, light_position);
            GL.Enable(EnableCap.Light0);

            // parameters for the camera
            Phi = 0.6f; Theta = 0.6f; Dist = 8.0f;
        }

        private void RenderingSettings(object sender, PaintEventArgs e)
        {
            // Modelview matrix
            GL.MatrixMode(MatrixMode.Modelview);
            GL.LoadIdentity();
            Matrix4 matLook = Matrix4.LookAt((float)(Dist * Math.Cos(Theta) * Math.Cos(Phi)), (float)(Dist * Math.Sin(Phi) * Math.Cos(Theta)), (float)(Dist * Math.Sin(Theta)), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
            GL.LoadMatrix(ref matLook);

            // perspective projection
            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();

            float aspect = 1;
            if((float)glControl.Height != 0 && (float)glControl.Width != 0) aspect = (float)glControl.Width / (float)glControl.Height;
            Matrix4  matPers = Matrix4.CreatePerspectiveFieldOfView(0.785f, aspect, 0.1f, 10.5f);
            
            GL.LoadMatrix(ref matPers);

            GL.Clear(ClearBufferMask.ColorBufferBit | ClearBufferMask.DepthBufferBit); 
            GL.Disable(EnableCap.Lighting);   
            GL.Enable(EnableCap.DepthTest);

            if (DrawCoordinateAxes) DrawAxes();  // draw coordinate axes

            GL.Enable(EnableCap.Lighting);
            GL.Enable(EnableCap.DepthTest);

            DrawTerrain(MainTerrain);     // draw terrain
            // GL.Disable(EnableCap.Lighting);

            if (DrawMesh)
            {
                //  GL.Enable(EnableCap.Lighting);
                // draw mesh
                GL.Translate(0.0f, 0.0f, 0.001f);  // posunutie droteneho modelu
                GL.PolygonMode(MaterialFace.Front, PolygonMode.Line);
                                
                GL.LineWidth(2.0f);  
                GL.Material(MaterialFace.Front, MaterialParameter.Specular, NetColor);
                GL.Material(MaterialFace.Front, MaterialParameter.Ambient, NetColor);

                GL.Begin(PrimitiveType.Triangles);
                for (int i = 0; i < MainTerrain.Indices.Count; i++)
                    GL.Vertex3(MainTerrain.Coordinates[MainTerrain.Indices[i]]);
                GL.End();
            }

            if (AddWaterLevel) DrawWaterLevel(); // draw water level

            // the buffers need to swapped, so the scene is drawn
            glControl.SwapBuffers();
        }

        private void RenderingAreaChanged(object sender, SizeChangedEventArgs e)
        {
            GL.Viewport(0, 0, glControl.Width, glControl.Height);
        }

        public Vector3 InverseProjection(Vector3 mouse, Matrix4 projection, Matrix4 view, System.Windows.Size viewport)
        {
            Vector4 vec;

            vec.X = 2.0f * mouse.X / (float)viewport.Width - 1;
            vec.Y = -(2.0f * mouse.Y / (float)viewport.Height - 1);
            vec.Z = mouse.Z;
            vec.W = 1.0f;

            Matrix4 viewInv = Matrix4.Invert(view);
            Matrix4 projInv = Matrix4.Invert(projection);

            Vector4.Transform(ref vec, ref projInv, out vec);
            Vector4.Transform(ref vec, ref viewInv, out vec);

            if (vec.W > 0.000001f || vec.W < -0.000001f)
            {
                vec.X /= vec.W;
                vec.Y /= vec.W;
                vec.Z /= vec.W;
            }

            return vec.Xyz;
        }

        //------------------------------------MOUSE MOVE, PRESS,... SETTINGS-----------------------------------------------------

        private void GLControl_MouseDown(object sender, System.Windows.Forms.MouseEventArgs e)
        {
            if (e.Button == MouseButtons.Right)
            {
                // change camera parameters
                IsRightDown = true;
                RightX = e.X;
                RightY = e.Y;
                oPhi = Phi;
                oTheta = Theta;
            }
            else if (e.Button == MouseButtons.Left)
            {
                // searching the control point beneath the mouse cursor 
                //when we are doing the inverse projection, what points lie in the ray which is casted from the point beneath the cursor
                // if there are any, we choose the closest one
                Vector3 start, end;
                int[] viewport = new int[4];
                Matrix4 modelMatrix, projMatrix;

                GL.GetFloat(GetPName.ModelviewMatrix, out modelMatrix);
                GL.GetFloat(GetPName.ProjectionMatrix, out projMatrix);
                GL.GetInteger(GetPName.Viewport, viewport);

                start = InverseProjection(new Vector3(e.X, e.Y, 0.0f), projMatrix, modelMatrix, new System.Windows.Size(viewport[2], viewport[3]));
                end = InverseProjection(new Vector3(e.X, e.Y, 1.0f), projMatrix, modelMatrix, new System.Windows.Size(viewport[2], viewport[3]));

                double se = Math.Sqrt(Vector3.Dot(start - end, start - end));

                for (int i = 0; i < MainTerrain.Coordinates.Count; i++)
                {
                    double sA = Math.Sqrt(Vector3.Dot(MainTerrain.Coordinates[i] - start, MainTerrain.Coordinates[i] - start));
                    double eA = Math.Sqrt(Vector3.Dot(MainTerrain.Coordinates[i] - end, MainTerrain.Coordinates[i] - end));

                    if (sA + eA > se - 0.001 && sA + eA < se + 0.001)
                    {
                        SelectedPoint = i;
                        IsLeftDown = true;

                        RightX = e.X;
                        RightY = e.Y;
                    }
                }
            }
        }

        private void GLControl_MouseMove(object sender, System.Windows.Forms.MouseEventArgs e)
        {
            if (IsRightDown) // camera rotation
            {
                IsRightDown = true;

                Phi = oPhi + (RightX - e.X) / 200.0f;
                Theta = oTheta + (e.Y - RightY) / 200.0f;
                glControl.Invalidate();  // redraw the scene
            }
            else if (IsLeftDown) // moving vertex
            {
                IsLeftDown = true;

                float Scaling = 0.003f;
                Vector3 staryBod = MainTerrain.Coordinates[SelectedPoint];
                MainTerrain.Coordinates[SelectedPoint] = new Vector3(staryBod.X, staryBod.Y, staryBod.Z + Convert.ToSingle(RightY - e.Y) * Scaling);
                int indexHeight = (int)Math.Floor((double)SelectedPoint / NumberOfPointsWidth);
                int indexWidth = SelectedPoint % NumberOfPointsWidth;
                MainTerrain.HeightmapAsArray[indexHeight, indexWidth] = MainTerrain.Coordinates[SelectedPoint].Z;
                MainTerrain.OldHeightmapAsArray[indexHeight, indexWidth] = MainTerrain.Coordinates[SelectedPoint].Z;


                RightY = e.Y;
                RightX = e.X;
                glControl.Invalidate(); // redraw the scene
            }
        }

        private void GLControl_MouseUp(object sender, System.Windows.Forms.MouseEventArgs e)
        {
            if (e.Button == MouseButtons.Right) IsRightDown = false;
            if (e.Button == MouseButtons.Left) IsLeftDown = false;
        }

        private void GLControl_MouseWheel(object sender, System.Windows.Forms.MouseEventArgs e)
        {
            Dist -= (double)e.Delta * 0.001; // zooming                        
            glControl.Invalidate(); // redraw the scene            
        }


        //------------------------------------------KEYBOARD SETTINGS------------------------------------------------------------

        private void Window_KeyDown(object sender, System.Windows.Input.KeyEventArgs e)
        {
            if (e.Key == Key.X) ChangeViewToYZ();
            if (e.Key == Key.Y) ChangeViewToXZ();
            if (e.Key == Key.Z) ChangeViewToXY();
        }

        private void ChangeViewToYZ()
        {
            if (IsX)
            {
                IsX = false;
                Phi = prevPhi;
                Theta = prevTheta;
                Dist = prevDist;
            }
            else
            {
                IsX = true;
                IsY = false;
                IsZ = false;
                prevPhi = Phi;
                prevTheta = Theta;
                prevDist = Dist;

                Phi = 0;
                Theta = 0.01;
                Dist = 8.0;
            }
            glControl.Invalidate();
        }

        private void ChangeViewToXZ()
        {
            if (IsY)
            {
                IsY = false;
                Phi = prevPhi;
                Theta = prevTheta;
                Dist = prevDist;
            }
            else
            {
                IsY = true;
                IsX = false;
                IsZ = false;
                prevPhi = Phi;
                prevTheta = Theta;
                prevDist = Dist;
                Phi = 1.57;
                Theta = 0.01;
                Dist = 8.0;
            }
            glControl.Invalidate();
        }

        private void ChangeViewToXY()
        {
            if (IsZ)
            {
                IsZ = false;
                Phi = prevPhi;
                Theta = prevTheta;
                Dist = prevDist;
            }
            else
            {
                IsZ = true;
                IsY = false;
                IsX = false;
                prevPhi = Phi;
                prevTheta = Theta;
                prevDist = Dist;
                Phi = 0;
                Theta = 1.57;
                Dist = 8.0;
            }
            glControl.Invalidate();
        }

        private void SaveWithShortcut(object sender, ExecutedRoutedEventArgs e)
        {
            SaveFile();
        }

        //----------------------------------------USER INTERFACE CONTROL---------------------------------------------------------

        //-----CHANGE OF VIEW-----
        private void ClickChangeToYZ(object sender, RoutedEventArgs e)
        {
            ChangeViewToYZ();
        }

        private void ClickChangeToXZ(object sender, RoutedEventArgs e)
        {
            ChangeViewToXZ();
        }

        private void ClickChangeToXY(object sender, RoutedEventArgs e)
        {
            ChangeViewToXY();
        }

        //--------SCALING---------
        private void XScaleSliderChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            if (MainTerrain != null)
            {
                MainTerrain.Scale.X = (float)XScaleSlider.Value;
                MainTerrain.RecomputeCoordinates();
            }
            glControl.Invalidate();
        }

        private void YScaleSliderChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            if (MainTerrain != null)
            {
                MainTerrain.Scale.Y = (float)YScaleSlider.Value;
                MainTerrain.RecomputeCoordinates();
            }
            glControl.Invalidate();
        }

        private void ZScaleSliderChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {

            if (MainTerrain != null)
            {
                float zScaleSliderValue = (float)ZScaleSlider.Value;
                MainTerrain.Scale.Z = zScaleSliderValue;
                // tieto veci su na vodu na jeho slider, sprav to pomocou min a max height, lebo takto je to zlozite
                float actualScale = zScaleSliderValue / oldZscale;
                MainTerrain.RecomputeCoordinates();
                if (zScaleSliderValue == 0) oldWaterLevelValue = (float)waterLevelSlider.Value;

                if (zScaleSliderValue >= 0)
                {
                    waterLevelSlider.Minimum = 0;
                    waterLevelSlider.Maximum = zScaleSliderValue;
                }
                else
                {
                    waterLevelSlider.Minimum = zScaleSliderValue;
                    waterLevelSlider.Maximum = 0;
                }

                if (Math.Sign(zScaleSliderValue) != Math.Sign(oldWaterLevelValue)) oldWaterLevelValue = -oldWaterLevelValue;
                waterLevelSlider.Value = oldZscale == 0 ? oldWaterLevelValue : actualScale * (float)waterLevelSlider.Value;
                oldZscale = zScaleSliderValue;
            }

            glControl.Invalidate();
        }

        //-------ITERATIONS-------
        private void IterationMinus(object sender, RoutedEventArgs e)
        {
            WasChangeOfIteration = true;
            if (NumberOfIterations > 0)
            {
                NumberOfPointsWidth = (int)Math.Ceiling((double)NumberOfPointsWidth / 2);
                NumberOfPointsHeight = (int)Math.Ceiling((double)NumberOfPointsHeight / 2);
                NumberOfIterations--;

                sliderWidthHeightmap.Value = NumberOfPointsWidth;
                sliderHeightHeightmap.Value = NumberOfPointsHeight;

                terrainInfoTextBlock.Text = "šírka terénu je " + NumberOfPointsWidth + " bodov a dĺžka je " + NumberOfPointsHeight + " bodov";
                terrainInfoTextBlock2.Text = "počet trojuholníkov siete je " + (2 * (NumberOfPointsWidth - 1) * (NumberOfPointsHeight - 1));
                iterationLabel.Content = Convert.ToString(NumberOfIterations);

                MainTerrain = new Terrain(MainTerrain.ReduceSampling(), NumberOfIterations, hurstSlider.Value, Scale);
                glControl.Invalidate(); // redraw the scene
            }
        }

        private void IterationPlus(object sender, RoutedEventArgs e)
        {
            WasChangeOfIteration = true;
            NumberOfPointsWidth = 2 * NumberOfPointsWidth - 1;
            NumberOfPointsHeight = 2 * NumberOfPointsHeight - 1;
            NumberOfIterations++;

            sliderWidthHeightmap.Value = NumberOfPointsWidth;
            sliderHeightHeightmap.Value = NumberOfPointsHeight;

            terrainInfoTextBlock.Text = "šírka terénu je " + NumberOfPointsWidth + " bodov a dĺžka je " + NumberOfPointsHeight + " bodov";
            terrainInfoTextBlock2.Text = "počet trojuholníkov siete je " + (2 * (NumberOfPointsWidth - 1) * (NumberOfPointsHeight - 1));
            iterationLabel.Content = Convert.ToString(NumberOfIterations);

            MainTerrain.Iteration = NumberOfIterations;
            MainTerrain = MainTerrain.DiamondSquareAlgorithm();
            glControl.Invalidate(); // redraw the scene
        }

        //-------WATER LEVEL-------
        private void WaterLevelSliderChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            glControl.Invalidate(); // redraw the scene
        }

        private void MouseDownOverWaterLevelPointsPicker(object sender, MouseButtonEventArgs e)
        {
            if (e.LeftButton == MouseButtonState.Pressed)
            {
                System.Windows.Point p = e.GetPosition(waterLevelPointsPicker);
                IndexOfMovedPointOfWaterLevel = -1;
                for (int i = 0; i < waterLevelPoints.Count; i++)
                    if (p.X >= waterLevelPoints[i].X - 5 && p.X <= waterLevelPoints[i].X + 5 && p.Y >= waterLevelPoints[i].Y - 5 && p.Y <= waterLevelPoints[i].Y + 5)
                        IndexOfMovedPointOfWaterLevel = i;

                if (IndexOfMovedPointOfWaterLevel != -1) return;

                waterLevelPoints.Add(p);
                if (waterLevelPoints.Count() > 3)
                {
                    if (IsIntersectionOfEdgesOfWaterLevel(waterLevelPoints.Count() - 1, waterLevelPoints))
                    {
                        System.Windows.MessageBox.Show(this, "Pri vytváraní vodenej hladiny sa prekrižujú hrany!", "Voľba bodov vodnej hladiny", MessageBoxButton.OK, MessageBoxImage.Error);
                        waterLevelPoints.RemoveAt(waterLevelPoints.Count() - 1);
                    }
                }
                DrawPolygonOverWaterLevelPointsPicker(waterLevelPointsPicker, waterLevelPoints);
            }

            if (e.RightButton == MouseButtonState.Pressed && waterLevelPoints.Count() > 3)
            {
                System.Windows.Point p = e.GetPosition(waterLevelPointsPicker);
                int nazvyInak = -1;
                for (int i = 0; i < waterLevelPoints.Count; i++)
                    if (p.X >= waterLevelPoints[i].X - 5 && p.X <= waterLevelPoints[i].X + 5 && p.Y >= waterLevelPoints[i].Y - 5 && p.Y <= waterLevelPoints[i].Y + 5)
                        nazvyInak = i;
                if (nazvyInak != -1) waterLevelPoints.RemoveAt(nazvyInak);
                DrawPolygonOverWaterLevelPointsPicker(waterLevelPointsPicker, waterLevelPoints);
            }
        }

        private void MouseMoveOverWaterLevelPointsPicker(object sender, System.Windows.Input.MouseEventArgs e)
        {
            if (e.LeftButton == MouseButtonState.Pressed && waterLevelPoints.Count() > 2)
            {
                System.Windows.Point p = e.GetPosition(waterLevelPointsPicker);
                if (p.X < 0) p.X = 0;
                if (p.X > waterLevelPointsPicker.Width) p.X = waterLevelPointsPicker.Width;
                if (p.Y < 0) p.Y = 0;
                if (p.Y > waterLevelPointsPicker.Height) p.Y = waterLevelPointsPicker.Height;
                if (IndexOfMovedPointOfWaterLevel != -1)
                {
                    waterLevelPoints[IndexOfMovedPointOfWaterLevel] = new System.Windows.Point(p.X, p.Y);
                    if (IsIntersectionOfEdgesOfWaterLevel(IndexOfMovedPointOfWaterLevel, waterLevelPoints))
                    {
                        System.Windows.MessageBox.Show(this, "Pri vytváraní vodenej hladiny sa prekrižujú hrany!", "Voľba bodov vodnej hladiny", MessageBoxButton.OK, MessageBoxImage.Error);
                        waterLevelPoints.RemoveAt(IndexOfMovedPointOfWaterLevel);
                    }
                }

                DrawPolygonOverWaterLevelPointsPicker(waterLevelPointsPicker, waterLevelPoints);
            }
        }

        private void WaterLevelClick(object sender, RoutedEventArgs e)
        {
            AddWaterLevel = !AddWaterLevel;
            waterLevelButton.Content = AddWaterLevel ? "Odober vodnú hladinu" : "Pridaj vodnú hladinu";
            waterLevelSlider.IsEnabled = AddWaterLevel;
            waterLevelPointsPicker.Children.Clear();
            waterLevelPointsPicker.IsEnabled = AddWaterLevel;
            waterLevelCutButton.IsEnabled = AddWaterLevel;
            waterLevelPoints = new List<System.Windows.Point>();
            IndexOfMovedPointOfWaterLevel = -1;
            glControl.Invalidate(); // redraw the scene
        }

        private void waterLevelCutButtonClick(object sender, RoutedEventArgs e)
        {
            ScanLine(MainTerrain);
            glControl.Invalidate(); // redraw the scene
        }

        //-----OTHER PARAMETERS----
        private void noiseSliderChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            if (MainTerrain != null)
            {
                float[,] heightmapAsArray = MainTerrain.HeightmapAsArray;
                if (PerlinNoiseActive) {
                    PerlinNoise simplexNoise = new PerlinNoise(512, 512, (float)noiseScaleSlider.Value, (int)noiseSeedSlider.Value); 
                    heightmapAsArray = ReduceHeightmap(ReduceHeightmap(simplexNoise.HeightmapAsArray));
                }

                if (SimplexNoiseActive) {
                    SimplexNoise simplexNoise = new SimplexNoise(512, 512, (float)noiseScaleSlider.Value, (int)noiseSeedSlider.Value); 
                    heightmapAsArray = ReduceHeightmap(ReduceHeightmap(simplexNoise.HeightmapAsArray));
                }
                
                MainTerrain = new Terrain(heightmapAsArray, MainTerrain.Iteration - 2, hurstSlider.Value, Scale);
                glControl.Invalidate(); // redraw the scene
            }            
        }

        private void hurstSliderChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            if (MainTerrain != null)
            {
                MainTerrain.Hurst = (float)hurstSlider.Value; 
                glControl.Invalidate(); // redraw the scene
            }
        }

        private void ClickFunctionFromCombobox(object sender, SelectionChangedEventArgs e)
        {
            if (MainTerrain != null)
            {
                PerlinNoiseActive = false;
                SimplexNoiseActive = false;
                MathFunctionActive = true;
                InitializeTerrainParameters(GetHeightmapFromFunction(functionInCombobox.SelectedIndex), ((ComboBoxItem)functionInCombobox.SelectedItem).Content.ToString());
                glControl.Invalidate(); // redraw the scene
            }
        }

        

        private void mySlider_DragStarted(object sender, System.Windows.Controls.Primitives.DragStartedEventArgs e)
        {
            sizeSliderDrag = true;
        }

        private void mySlider_DragCompleted(object sender, System.Windows.Controls.Primitives.DragCompletedEventArgs e)
        {
            sizeSliderDrag = false;
        }


        private void WidthHeightSliderChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            if (sizeSliderDrag)
            {
                if (MainTerrain != null)
                {
                    float[,] newHeightmap = new float[(int)sliderHeightHeightmap.Value, (int)sliderWidthHeightmap.Value];
                    if (MathFunctionActive == false & PerlinNoiseActive == false && SimplexNoiseActive == false)
                    {
                        newHeightmap = ChangeHeightmapSize();
                    }
                    else {
                        if (!WasChangeOfIteration)
                        {
                            if (MathFunctionActive)
                            {
                                newHeightmap = GetHeightmapFromFunction(functionInCombobox.SelectedIndex);
                            }
                            if (PerlinNoiseActive)
                            {
                                PerlinNoise perlinNoise = new PerlinNoise((int)sliderWidthHeightmap.Value * 4, (int)sliderHeightHeightmap.Value * 4, (float)noiseScaleSlider.Value, (int)noiseSeedSlider.Value);
                                newHeightmap = ReduceHeightmap(ReduceHeightmap(perlinNoise.HeightmapAsArray));
                            }

                            if (SimplexNoiseActive)
                            {
                                SimplexNoise simplexNoise = new SimplexNoise((int)sliderWidthHeightmap.Value * 4, (int)sliderHeightHeightmap.Value * 4, (float)noiseScaleSlider.Value, (int)noiseSeedSlider.Value);
                                newHeightmap = ReduceHeightmap(ReduceHeightmap(simplexNoise.HeightmapAsArray));
                            }
                        }
                        else
                        {
                            newHeightmap = ChangeHeightmapSize();
                        }
                    }

                    NumberOfPointsWidth = newHeightmap.GetLength(1);
                    NumberOfPointsHeight = newHeightmap.GetLength(0);

                    NumberOfIterations = NumberOfPointsWidth >= NumberOfPointsHeight ? (int)Math.Log(NumberOfPointsWidth, 2) : (int)Math.Log(NumberOfPointsHeight, 2);
                    terrainInfoTextBlock.Text = "šírka terénu je " + NumberOfPointsWidth + " bodov a dĺžka je " + NumberOfPointsHeight + " bodov";
                    terrainInfoTextBlock2.Text = "počet trojuholníkov siete je " + (2 * (NumberOfPointsWidth - 1) * (NumberOfPointsHeight - 1));
                    iterationLabel.Content = Convert.ToString(NumberOfIterations);

                    MainTerrain = new Terrain(newHeightmap, NumberOfIterations, hurstSlider.Value, Scale);                    
                    glControl.Invalidate(); // redraw the scene

                }
            }
        }

        private float [,] ChangeHeightmapSize()
        {
            float[,] newHeightmap = new float[(int)sliderHeightHeightmap.Value, (int)sliderWidthHeightmap.Value];
            if (newHeightmap.GetLength(1) > MainTerrain.HeightmapAsArray.GetLength(1) || newHeightmap.GetLength(0) > MainTerrain.HeightmapAsArray.GetLength(0))
            {
                if (newHeightmap.GetLength(1) > MainTerrain.HeightmapAsArray.GetLength(1))
                {
                    for (int i = 0; i < newHeightmap.GetLength(0); i++)
                    {
                        for (int j = 0; j < newHeightmap.GetLength(1); j++)
                        {
                            if (j > MainTerrain.HeightmapAsArray.GetLength(1) - 1)
                            {
                                newHeightmap[i, j] = MainTerrain.HeightmapAsArray[i, MainTerrain.HeightmapAsArray.GetLength(1) - 1];
                            }
                            else
                            {
                                newHeightmap[i, j] = MainTerrain.HeightmapAsArray[i, j];
                            }
                        }
                    }
                }

                if (newHeightmap.GetLength(0) > MainTerrain.HeightmapAsArray.GetLength(0))
                {
                    for (int i = 0; i < newHeightmap.GetLength(0); i++)
                    {
                        for (int j = 0; j < newHeightmap.GetLength(1); j++)
                        {
                            if (i > MainTerrain.HeightmapAsArray.GetLength(0) - 1)
                            {
                                newHeightmap[i, j] = MainTerrain.HeightmapAsArray[MainTerrain.HeightmapAsArray.GetLength(0) - 1, j];

                            }
                            else
                            {
                                newHeightmap[i, j] = MainTerrain.HeightmapAsArray[i, j];
                            }

                        }
                    }
                }
            }
            else
            {
                for (int i = 0; i < newHeightmap.GetLength(0); i++)
                {
                    for (int j = 0; j < newHeightmap.GetLength(1); j++)
                    {
                        newHeightmap[i, j] = MainTerrain.HeightmapAsArray[i, j];
                    }
                }
            }



            return newHeightmap;
        }


        private void FunctionSliderChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            if (MainTerrain != null)
            {
                InitializeTerrainParameters(GetHeightmapFromFunction(functionInCombobox.SelectedIndex), functionInCombobox.Text);
                glControl.Invalidate(); // redraw the scene
            }             
        }

        //-----MENU PARAMETERS-----
        private void ClickNew(object sender, RoutedEventArgs e)
        {
            StartWindow startwindow = new StartWindow();
            startwindow.Show();
            Close();
        }

        private void ClickOpenFile(object sender, RoutedEventArgs e)
        {
            OpenFile();
        }

        private void ClickSaveFile(object sender, RoutedEventArgs e)
        {
            SaveFile();
        }

        private void ClickExit(object sender, RoutedEventArgs e)
        {
            Close();
        }

        private void ClickShowAxes(object sender, RoutedEventArgs e)
        {
            DrawCoordinateAxes = !DrawCoordinateAxes;
            glControl.Invalidate(); // redraw the scene
        }

        private void ClickShowMesh(object sender, RoutedEventArgs e)
        {
            DrawMesh = !DrawMesh;
            glControl.Invalidate(); // redraw the scene

        }

        private void ClickManual(object sender, RoutedEventArgs e)
        {
            string information = "V pravej časti aplikácie sa nastavujú všetky parametre, ktoré môžeme ovládať pri teréne. " 
                + "Vo vrchnej časti sa nachádzajú všeobecné parametre terénu ako škálovanie a pridávanie/uberanie bodov terénu a nachádza sa tu tlačidlo na vyhladenie terénu. "
                + "Treba si uvedomiť, že škálovanie v smere osí x a y sa na uložených výškových mapách neprejaví, kedže tie majú vždy pravidelné rozostupy. Škálovanie v týchto smeroch je užitočné pre exportovanie 3D modelu. \n"
                + "Následne sú časti parametrov rozdelené do celkov. Máme tu parametre na ovládanie šumu a parametre pre zmenu matematických funkcií. \n"
                + "Dôležitým parametrom v tejto časti je zmena šírky a dĺžky terénu. Pri pridávaní, v ktoromkoľvek smere sa kopírujú okrajové hodnoty terénu. Ak terén zmenšujeme tak sa iba uberajú posledné riadky/stĺpce výškovej mapy."
                + "Tento parameter je zaujímavý pri funkciách a šumoch, len pri nich sa najskôr správa inak. Najskôr ním definujeme veľkosť šumu, prípadne počet vzorkovaných bodov funkcie a keď už sa zmení počet iterácii terénu (teda už s ním nejako pracujeme)," 
                + "tak sa znova správa, akoby sme mali bežný terén. \n"
                + "V spodnej časti sa nachádza rozhranie na pridávanie vodnej hladiny. Je tam plocha, do ktorej sa zadáva obrys vodnej hladiny, pri ktorom ale nemôžeme vytvárať prekrývajúce sa hrany. "
                + "Po vytvorení obrysu vodnej hladiny jej vieme nastaviť výšku a následne tlačidlom na orezanie terénu môžeme terén upraviť pomocou zvolenej hladiny. \n" 
                + "\n"
                + "Na vrchnej lište menu si vieme zvoliť nový typ terénu, šumu alebo matematickej funkcie. V časti Zobraziť môžeme zmeniť pohľad na terén a povoliť vykresľovanie pletiva a súradnicových osí. \n"
                + "Časť súbor nám ponúka otvorenie úvodu aplikácie a uloženie a načítanie terénu.";
            System.Windows.MessageBox.Show(information, "Informácie", MessageBoxButton.OK, MessageBoxImage.Information);
        }

        private void ClickTerrainFromMenu(object sender, RoutedEventArgs e)
        {
            PerlinNoiseActive = false;
            SimplexNoiseActive = false;
            MathFunctionActive = false;
            System.Windows.Controls.MenuItem selectedTerrainItem = sender as System.Windows.Controls.MenuItem;
            string terrainType = selectedTerrainItem.Name;
            string applicationPath = Environment.CurrentDirectory;
            string cesta = applicationPath.Remove(applicationPath.Length - 9) + "Heightmaps/" + terrainType + ".png"; //leave actual path and open file Heightmaps/
            Bitmap heightmapAsImage = (Bitmap)System.Drawing.Image.FromFile(cesta, true);

            InitializeTerrainParameters(GetHeightmapFromImage(heightmapAsImage), terrainType);
            glControl.Invalidate(); // redraw the scene                                    
        }

        private void ClickPerlinFromMenu(object sender, RoutedEventArgs e)
        {
            PerlinNoiseActive = true;
            SimplexNoiseActive = false;
            MathFunctionActive = false;
            PerlinNoise perlinNoise = new PerlinNoise((int)sliderWidthHeightmap.Value * 4, (int)sliderHeightHeightmap.Value * 4, (float)noiseScaleSlider.Value, (int)noiseSeedSlider.Value);
            float[,] heightmapAsArray = ReduceHeightmap(ReduceHeightmap(perlinNoise.HeightmapAsArray));
            string terrainType = "perlinov šum";
            InitializeTerrainParameters(heightmapAsArray, terrainType);
            glControl.Invalidate(); // redraw the scene
        }

        private void ClickSimplexFromMenu(object sender, RoutedEventArgs e)
        {
            PerlinNoiseActive = false;
            SimplexNoiseActive = true;
            MathFunctionActive = false;
            SimplexNoise simplexNoise = new SimplexNoise((int)sliderWidthHeightmap.Value * 4, (int)sliderHeightHeightmap.Value * 4, (float)noiseScaleSlider.Value, (int)noiseSeedSlider.Value);
            float[,] heightmapAsArray = ReduceHeightmap(ReduceHeightmap(simplexNoise.HeightmapAsArray));
            string terrainType = "simplexný šum";
            InitializeTerrainParameters(heightmapAsArray, terrainType);
            glControl.Invalidate(); // redraw the scene
        }

        private void ClickFunctionFromMenu(object sender, RoutedEventArgs e)
        {
            PerlinNoiseActive = false;
            SimplexNoiseActive = false;
            MathFunctionActive = true;
            System.Windows.Controls.MenuItem selectedFunctionItem = sender as System.Windows.Controls.MenuItem;
            functionInCombobox.SelectedItem = functionInCombobox.Items[functionInMenu.Items.IndexOf(selectedFunctionItem)];
            InitializeTerrainParameters(GetHeightmapFromFunction(functionInMenu.Items.IndexOf(selectedFunctionItem)), selectedFunctionItem.Header.ToString());
            glControl.Invalidate(); // redraw the scene
        }

        //---------------------------------------------FUNCTIONS----------------------------------------------------------

        //------OPEN AND SAVE------
        private void OpenFile()
        {
            string filePath = null;
            Microsoft.Win32.OpenFileDialog openFileDialog = new Microsoft.Win32.OpenFileDialog();
            openFileDialog.FileName = "";
            openFileDialog.Filter = "Súbory bitových máp (*.bmp)|*.bmp" +
                "|JPEG (*.jpg;*.jpeg;*.jpe;*.jfif)|*.jpg;*.jpeg;*.jpe;*.jfif" +
                "|GIF (*.gif)|*.gif" +
                "|TIFF (*.tif;*.tiff)|*.tif;*.tiff" +
                "|PNG (*.png)|*.png" +
                "|ICO (*.ico)|*.ico" +
                "|Všetky súbory|*.*";
            openFileDialog.FilterIndex = 7;

            if (openFileDialog.ShowDialog() == true)
            {
                filePath = openFileDialog.FileName;
            }

            if (filePath == null)
                return;
            try
            {
                Bitmap image = (Bitmap)System.Drawing.Image.FromFile(filePath, true);
                float[,] mapa = GetHeightmapFromImage(image);
                InitializeTerrainParameters(mapa, "Vlastný obrázok");
                glControl.Invalidate();  // redraw the scene
            }
            catch (Exception a)
            {
                System.Windows.MessageBox.Show(this, "Nepodporovaný súbor!", "Otvorenie súboru", MessageBoxButton.OK, MessageBoxImage.Error);
                return;
            }
        }

        private void SaveFile()
        {
            int index = -1;
            string filePath = null;
            Microsoft.Win32.SaveFileDialog saveFileDialog = new Microsoft.Win32.SaveFileDialog();
            saveFileDialog.Filter = "Súbory bitových máp (*.bmp)|*.bmp" +
                "|JPEG (*.jpg;*.jpeg;*.jpe;*.jfif)|*.jpg;*.jpeg;*.jpe;*.jfif" +
                "|GIF (*.gif)|*.gif" +
                "|TIFF (*.tif;*.tiff)|*.tif;*.tiff" +
                "|PNG (*.png)|*.png" +
                "|OBJ (*.obj)|*.obj";
            saveFileDialog.FilterIndex = 1;
            if (saveFileDialog.ShowDialog() == true)
            {
                filePath = saveFileDialog.FileName;
                index = saveFileDialog.FilterIndex;
            }

            if (filePath == null || index == -1)
                return;

            try
            {
                System.Drawing.Imaging.ImageFormat imageformat = null;

                switch (index)
                {
                    case 1:
                        imageformat = System.Drawing.Imaging.ImageFormat.Bmp; break;
                    case 2:
                        imageformat = System.Drawing.Imaging.ImageFormat.Jpeg; break;
                    case 3:
                        imageformat = System.Drawing.Imaging.ImageFormat.Gif; break;
                    case 4:
                        imageformat = System.Drawing.Imaging.ImageFormat.Tiff; break;
                    case 5:
                        imageformat = System.Drawing.Imaging.ImageFormat.Png; break;
                }

                if (index > 0 && index < 6)
                {
                    SaveArrayAsHeightmap(MainTerrain.HeightmapAsArray, filePath, imageformat);
                }
                else
                {
                    SaveArrayAsOBJ(MainTerrain, filePath);
                }
            }
            catch (Exception a)
            {
                System.Windows.MessageBox.Show(this, "Nastala chyba " + a, "Chybová hláška", MessageBoxButton.OK, MessageBoxImage.Error);
                return;
            }
        }

        private void SaveArrayAsOBJ(Terrain terrain, string path)
        {
            string vysledok = "";
            for (int i = 0; i < terrain.Coordinates.Count; i++)
            {
                vysledok += "v " + terrain.Coordinates[i].X + " " + (terrain.Coordinates[i].Z) + " " + terrain.Coordinates[i].Y + "\n";
            }

            for (int i = 0; i < terrain.Indices.Count; i += 3)
            {
                vysledok += "f " + (terrain.Indices[i] + 1) + " " + (terrain.Indices[i + 1] + 1) + " " + (terrain.Indices[i + 2] + 1) + "\n";
            }

            File.WriteAllText(path, vysledok);
        }

        private void SaveArrayAsHeightmap(float[,] heightmapAsArray, string path, System.Drawing.Imaging.ImageFormat imageFormat)
        {
            int width = heightmapAsArray.GetLength(1);
            int height = heightmapAsArray.GetLength(0);
            Bitmap bmp = new Bitmap(width, height);
            Vector2 minAndMaxHeight = GetMinAndMaxHeightFromArray(heightmapAsArray);  // X is min value and Y is max value
            float interval = minAndMaxHeight.Y - minAndMaxHeight.X;
            if (interval < 1) interval = 1;
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    int color = (int)Math.Round((MainTerrain.HeightmapAsArray[i, j] - minAndMaxHeight.X) / interval * 255);  // normalization data to interval <0,1> and then creating color
                    bmp.SetPixel(j, i, System.Drawing.Color.FromArgb(color, color, color));
                }
            }

            var stream = File.OpenWrite(path);
            bmp.Save(stream, imageFormat);
            stream.Close();
            bmp.Dispose();
        }

        //------WATER LEVEL--------
        private void DrawWaterLevel()
        {            
            double waterLevelHeight = waterLevelSlider.Value + 0.005f;

            float[] diffuse = { 0.9f, 0.9f, 0.9f, 1.0f };
            float[] specular = { 0.1f, 0.1f, 0.1f, 0.5f };

            GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Diffuse, diffuse);
            GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Specular, specular);
            GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Shininess, 0.1f);
            GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Ambient, WaterColor);
            GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill);

            GL.Begin(PrimitiveType.Polygon);
            for (int i = 0; i < waterLevelPoints.Count(); i++)
            {
                Vector3 bod = MainTerrain.Coordinates[GetWaterLevelCoordinateIndex(waterLevelPoints[i])];
                GL.Normal3(0.0f, 0.0f, 1.0f);
                GL.Vertex3(bod.X, bod.Y, waterLevelHeight);
            }

            GL.End();
        }

        private void DrawPolygonOverWaterLevelPointsPicker(Canvas canvas, List<System.Windows.Point> listOfPoints)
        {
            canvas.Children.Clear();
            if (listOfPoints.Count() > 2)
            {
                Polygon myPolygon = new Polygon();
                myPolygon.Stroke = System.Windows.Media.Brushes.Black;
                myPolygon.Fill = new SolidColorBrush(System.Windows.Media.Color.FromScRgb(1.0f, WaterColor[0], WaterColor[1], WaterColor[2]));
                //myPolygon.FillRule = FillRule.Nonzero;
                myPolygon.Points = new PointCollection(listOfPoints);
                canvas.Children.Add(myPolygon);
            }

            if (listOfPoints.Count() > 2)
            {
                System.Windows.Vector vector = listOfPoints[listOfPoints.Count() - 2] - listOfPoints[listOfPoints.Count() - 1];
                vector.Normalize();
                for (int k = 0; k < 2; k++)
                {
                    Line l = new Line();
                    l.Stroke = new SolidColorBrush(Colors.Black);
                    l.StrokeThickness = 1;
                    l.X1 = listOfPoints[listOfPoints.Count() - 1].X;
                    l.Y1 = listOfPoints[listOfPoints.Count() - 1].Y;
                    l.X2 = listOfPoints[listOfPoints.Count() - 1].X + (vector.X * 0.866 - vector.Y * 0.5) * 10;
                    l.Y2 = listOfPoints[listOfPoints.Count() - 1].Y + (vector.X * 0.5 + vector.Y * 0.866) * 10;
                    if (k == 1)
                    {
                        l.X2 = listOfPoints[listOfPoints.Count() - 1].X + (vector.X * 0.866 + vector.Y * 0.5) * 10;
                        l.Y2 = listOfPoints[listOfPoints.Count() - 1].Y + (-vector.X * 0.5 + vector.Y * 0.866) * 10;
                    }
                    waterLevelPointsPicker.Children.Add(l);
                }
            }

            for (int i = 0; i < listOfPoints.Count(); i++)
            {
                Ellipse myEllipse = new Ellipse();
                myEllipse.Fill = new SolidColorBrush(Colors.Black);
                if (i == listOfPoints.Count() - 1) myEllipse.Fill = new SolidColorBrush(Colors.Red);
                myEllipse.Width = 4;
                myEllipse.Height = 4;
                Canvas.SetLeft(myEllipse, listOfPoints[i].X - 2);
                Canvas.SetTop(myEllipse, listOfPoints[i].Y - 2);
                waterLevelPointsPicker.Children.Add(myEllipse);
            }

            glControl.Invalidate(); // redraw the scene
        }

        private int GetIndexOnInterval(double percentualExpressionOfCoordinate, int dimensionLength)
        {
            int bottomIndex = 0;
            int topIndex = (dimensionLength - 1);

            while (bottomIndex != topIndex)
            {
                int i = (int)Math.Floor((double)(topIndex + bottomIndex) / 2);
                if (percentualExpressionOfCoordinate < (double)i / (dimensionLength - 1))
                {
                    topIndex = i;
                }
                else
                {
                    bottomIndex = i;
                }

                if ((topIndex - bottomIndex) == 1)
                {
                    if ((percentualExpressionOfCoordinate - (double)bottomIndex / (dimensionLength - 1)) < ((double)topIndex / (dimensionLength - 1) - percentualExpressionOfCoordinate))
                    {
                        topIndex = bottomIndex;
                    }
                    else
                    {
                        bottomIndex = topIndex;
                    }
                }
            }
            return bottomIndex;
        }

        private int GetWaterLevelCoordinateIndex(System.Windows.Point point)
        {
            int widthIndex = GetIndexOnInterval(point.X / waterLevelPointsPicker.Width, NumberOfPointsWidth);
            int heightIndex = GetIndexOnInterval(point.Y / waterLevelPointsPicker.Height, NumberOfPointsHeight);
            return heightIndex * NumberOfPointsWidth + widthIndex;
        }


        //------OTHERS FUNCTIONS-----


        private void DrawAxes()
        {
            GL.Begin(PrimitiveType.Lines);
            GL.Enable(EnableCap.Lighting);
            GL.Disable(EnableCap.DepthTest);

            GL.Color3(1.0f, 0.0f, 0.0f);
            GL.Vertex3(0.0f, 0.0f, 0.0f);
            GL.Color3(1.0f, 0.0f, 0.0f);
            GL.Vertex3(20.0f, 0.0f, 0.0f);

            GL.Color3(0.0f, 1.0f, 0.0f);
            GL.Vertex3(0.0f, 0.0f, 0.0f);
            GL.Color3(0.0f, 1.0f, 0.0f);
            GL.Vertex3(0.0f, 20.0f, 0.0f);

            GL.Color3(0.0f, 0.0f, 1.0f);
            GL.Vertex3(0.0f, 0.0f, 0.0f);
            GL.Color3(0.0f, 0.0f, 1.0f);
            GL.Vertex3(0.0f, 0.0f, 20.0f);
            GL.End();

        }

        private void DrawTerrain(Terrain terrain)
        {            
            GL.PolygonMode(MaterialFace.Front, PolygonMode.Fill); // enabble filling of shapes with color 

            // edit if you want something different
            float[] diffuse = { 0.9f, 0.9f, 0.9f, 1.0f };
            float[] specular = { 0.1f, 0.1f, 0.1f, 0.5f };

            GL.Material(MaterialFace.Front, MaterialParameter.Diffuse, diffuse);
            GL.Material(MaterialFace.Front, MaterialParameter.Specular, specular);
            GL.Material(MaterialFace.Front, MaterialParameter.Shininess, 0.1f);

            PrimitiveType netType = PrimitiveType.Triangles;

            // draw triangles 
            GL.Begin(netType);
            for (int i = 0; i < terrain.Indices.Count; i++)
            {
                GL.Material(MaterialFace.Front, MaterialParameter.Ambient, TerrainColor);
                GL.Normal3(0.0f, 0.0f, 1.0f); 
                GL.Vertex3(terrain.Coordinates[terrain.Indices[i]]);
            }
            GL.End();            
        }

        public float[,] ReduceHeightmap(float[,] heightmap)
        {
            int oldWidth = heightmap.GetLength(1);
            int oldHeight = heightmap.GetLength(0);

            int newWidth = (int)Math.Ceiling((float)oldWidth / 2);
            int newHeight = (int)Math.Ceiling((float)oldHeight / 2);
            if (newWidth < 2) newWidth = oldWidth;
            if (newHeight < 2) newHeight = oldHeight;
            float[,] newHeightmap = new float[newHeight, newWidth];

            for (int i = 0; i < oldHeight; i = i + 2)
                for (int j = 0; j < oldWidth; j = j + 2)
                {
                    newHeightmap[i / 2, j / 2] = heightmap[i, j];
                }
            return newHeightmap;
        }

        public float[,] Convolution(float[,] heightmap)
        {
            int[,] convolutionMatrix = new int[3, 3] { { 1, 2, 1 }, { 2, 4, 2 }, { 1, 2, 1 } };
            float[,] newHeightmap = (float[,])heightmap.Clone();
            float newValue;
            Vector2 minAndMaxValue = GetMinAndMaxHeightFromArray(heightmap);

            for (int i = 0; i < heightmap.GetLength(0) - 2; i++)
            {
                for (int j = 0; j < heightmap.GetLength(1) - 2; j++)
                {
                    newValue = (float)(((heightmap[i, j] * convolutionMatrix[0, 0] +
                                 (heightmap[i, j + 1] * convolutionMatrix[0, 1]) +
                                 (heightmap[i, j + 2] * convolutionMatrix[0, 2]) +
                                 (heightmap[i + 1, j] * convolutionMatrix[1, 0]) +
                                 (heightmap[i + 1, j + 1] * convolutionMatrix[1, 1]) +
                                 (heightmap[i + 1, j + 2] * convolutionMatrix[1, 2]) +
                                 (heightmap[i + 2, j] * convolutionMatrix[2, 0]) +
                                 (heightmap[i + 2, j + 1] * convolutionMatrix[2, 1]) +
                                 (heightmap[i + 2, j + 2] * convolutionMatrix[2, 2])) * 0.0625));

                    if (newValue < minAndMaxValue.X)
                    {
                        newValue = minAndMaxValue.X;
                    }
                    else if (newValue >= minAndMaxValue.Y)
                    {
                        newValue = minAndMaxValue.Y;
                    }

                    newHeightmap[i + 1, j + 1] = newValue;
                }
            }
            return newHeightmap;
        }

        private Vector2 GetMinAndMaxHeightFromBitmap(Bitmap image1)
        {
            int maxHeight = -1;
            int minHeight = 256;

            for (int i = 0; i < image1.Height; i++)
            {
                for (int j = 0; j < image1.Width; j++)
                {
                    if (image1.GetPixel(j, i).R >= maxHeight) maxHeight = image1.GetPixel(j, i).R;
                    if (image1.GetPixel(j, i).R <= minHeight) minHeight = image1.GetPixel(j, i).R;
                }
            }
            return new Vector2(minHeight, maxHeight);
        }

        private Vector2 GetMinAndMaxHeightFromArray(float[,] heightmapAsArray)
        {
            float maxHeight = float.MinValue;
            float minHeight = float.MaxValue;

            for (int i = 0; i < heightmapAsArray.GetLength(0); i++)
            {
                for (int j = 0; j < heightmapAsArray.GetLength(1); j++)
                {
                    if (heightmapAsArray[i, j] >= maxHeight) maxHeight = heightmapAsArray[i, j];
                    if (heightmapAsArray[i, j] <= minHeight) minHeight = heightmapAsArray[i, j];
                }
            }
            return new Vector2(minHeight, maxHeight);
        }

        private float[,] GetHeightmapFromFunction(int functionIndex)
        {
            int width = (int)sliderWidthHeightmap.Value;
            int height = (int)sliderHeightHeightmap.Value;
            float minX = (float)rangeSliderMinX.Value;
            float minY = (float)rangeSliderMinY.Value;
            float rangeX = (float)Math.Abs(rangeSliderMaxX.Value - minX);
            float rangeY = (float)Math.Abs(rangeSliderMaxY.Value - minY);
            float[,] heightMap = new float[height, width];
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                {
                    float x = minX + rangeX * i / (height - 1);
                    float y = minY + rangeY * j / (width - 1);

                    switch (functionIndex)
                    {
                        case 0:
                            heightMap[i, j] = (float)(Math.Sin(x + y) + 1.0f) / 2;
                            break;
                        case 1:
                            heightMap[i, j] = (float)(Math.Sin(Math.Cos(x - y)) + 1.0f) / 2;
                            break;
                        case 2:
                            heightMap[i, j] = (float)((Math.Cos(x) * Math.Sin(y)) + 1.0f) / 2;
                            break;
                        case 3:
                            heightMap[i, j] = (float)(Math.Cos(x) - Math.Sin(y) + 2.0f) / 4;
                            break;
                        case 4:
                            heightMap[i, j] = (float)(Math.Cos(x) * Math.Cos(x) * Math.Cos(y) * Math.Cos(y));
                            break;
                        case 5:
                            heightMap[i, j] = (float)Math.Pow(Math.E, -(x * x + y * y));
                            break;
                        case 6:
                            heightMap[i, j] = (float)((x * x + y * y) * Math.Pow(Math.E, 1 - x * x - y * y));
                            break;
                        case 7:
                            heightMap[i, j] = (float)(Math.Sin(x * x + y * y) + 1.0f) / 2;
                            break;
                        case 8:
                            heightMap[i, j] = (float)(Math.Sin(x * x + y * y) / (x * x + y * y) * 0.81 + 0.19);
                            break;
                    }
                }
            return heightMap;
        }

        private float[,] GetHeightmapFromImage(Bitmap image)
        {
            int width = image.Width;
            int height = image.Height;
            Vector2 minimumAndMaximumHeight = GetMinAndMaxHeightFromBitmap(image);  // first coordinate is minimum and second is maximum height
            float[,] heightmapAsArray = new float[height, width]; 
            float interval = minimumAndMaximumHeight.Y - minimumAndMaximumHeight.X;
            if (interval == 0) interval = 1;
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {                    
                    heightmapAsArray[i, j] = (float)(image.GetPixel(j, i).R - minimumAndMaximumHeight.X) / 255; 

                }
            }
            return heightmapAsArray;
        }

        private void sliderWidthHeightmap_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
           // Console.WriteLine("Mame tam");
        }
        private void drag(object sender, System.Windows.Controls.Primitives.DragStartedEventArgs e)
        {
            Console.WriteLine("Tahame slider");
        }


        private void TerrainSmoothingButtonClick(object sender, RoutedEventArgs e)
        {
            MainTerrain = new Terrain(Convolution(MainTerrain.HeightmapAsArray), MainTerrain.Iteration, MainTerrain.Hurst, MainTerrain.Scale);
            glControl.Invalidate();  // redraw the scene
        }
                
        private bool IsIntersectionOfEdgesOfWaterLevel(int index, List<System.Windows.Point> listOfPoints)
        {
            int prevIndex = listOfPoints.Count() - 1;
            int nextIndex = 0;
            if ((index + 1) < listOfPoints.Count()) nextIndex = index + 1;
            if ((index - 1) >= 0) prevIndex = index - 1;

            bool intersection = false;

            int stop = prevIndex - 1;
            if (stop < 0) { stop = listOfPoints.Count() - 1; }
            int k = nextIndex;
            while (k != stop)
            {
                int nextI = k + 1;
                if (nextI == listOfPoints.Count()) nextI = 0;
                if (doEdgesIntersect(listOfPoints[index], listOfPoints[prevIndex], listOfPoints[k], listOfPoints[nextI])) intersection = true;
                k = nextI;
            }

            k = nextIndex + 1;
            if (k == listOfPoints.Count()) k = 0;
            while (k != prevIndex)
            {
                int nextI = k + 1;
                if (nextI == listOfPoints.Count()) nextI = 0;
                if (doEdgesIntersect(listOfPoints[index], listOfPoints[nextIndex], listOfPoints[k], listOfPoints[nextI])) intersection = true;
                k = nextI;
            }
            return intersection;
        }

        static bool IsPointOnEdge(System.Windows.Point p, System.Windows.Point q, System.Windows.Point r)
        {
            if (q.X <= Math.Max(p.X, r.X) && q.X >= Math.Min(p.X, r.X) &&
                q.Y <= Math.Max(p.Y, r.Y) && q.Y >= Math.Min(p.Y, r.Y))
                return true;

            return false;
        }

        static int WaterLevelPolygonOrientation(System.Windows.Point p, System.Windows.Point q, System.Windows.Point r)
        {
            double val = (q.Y - p.Y) * (r.X - q.X) - (q.X - p.X) * (r.Y - q.Y);
            if (val == 0) return 0; // collinear
            return (val > 0) ? 1 : 2; // clock or counterclock wise
        }

        static bool doEdgesIntersect(System.Windows.Point p1, System.Windows.Point q1, System.Windows.Point p2, System.Windows.Point q2)
        {
            // Find the four orientations needed for general and special cases
            int o1 = WaterLevelPolygonOrientation(p1, q1, p2);
            int o2 = WaterLevelPolygonOrientation(p1, q1, q2);
            int o3 = WaterLevelPolygonOrientation(p2, q2, p1);
            int o4 = WaterLevelPolygonOrientation(p2, q2, q1);

            if (o1 != o2 && o3 != o4)
                return true;
                        
            if (o1 == 0 && IsPointOnEdge(p1, p2, q1)) return true;            
            if (o2 == 0 && IsPointOnEdge(p1, q2, q1)) return true;
            if (o3 == 0 && IsPointOnEdge(p2, p1, q2)) return true;
            if (o4 == 0 && IsPointOnEdge(p2, q1, q2)) return true;
            return false; 
        }

        private List<System.Windows.Point> GetMinMaxBoxOfWaterLevel(List<System.Windows.Point> waterLevelPoints)
        {
            List<System.Windows.Point> minMaxBox = new List<System.Windows.Point>();
            float minX = float.MaxValue;
            float minY = float.MaxValue;
            float maxX = float.MinValue;
            float maxY = float.MinValue;

            foreach (System.Windows.Point p in waterLevelPoints)
            {
                if (p.X < minX) minX = (float)p.X;
                if (p.X > maxX) maxX = (float)p.X;
                if (p.Y < minY) minY = (float)p.Y;
                if (p.Y > maxY) maxY = (float)p.Y;
            }

            minMaxBox.Add(new System.Windows.Point(minX, minY));
            minMaxBox.Add(new System.Windows.Point(maxX, maxY));
            return minMaxBox;
        }

        private void SetHeigntOnInterval(Vector3 startPoint, Vector3 endPoint, Terrain terrain)
        {
            List<System.Windows.Point> minMaxBox = GetMinMaxBoxOfWaterLevel(waterLevelPoints);
            int startIndexOfBox = GetWaterLevelCoordinateIndex(minMaxBox[0]);
            int endIndexOfBox = GetWaterLevelCoordinateIndex(minMaxBox[1]);
            float waterLevelHeight = (float)waterLevelSlider.Value;
            int minMaxBoxWidthPoints = (endIndexOfBox % NumberOfPointsWidth) - (startIndexOfBox % NumberOfPointsWidth);  // the number of points per width in the box

            for (int firstElementOfRow = startIndexOfBox; firstElementOfRow <= endIndexOfBox; firstElementOfRow += NumberOfPointsWidth)
            {
                if (terrain.Coordinates[firstElementOfRow].X < startPoint.X) continue;
                if (terrain.Coordinates[firstElementOfRow].X > startPoint.X) break;
                for (int columnIndexInRow = firstElementOfRow; columnIndexInRow <= firstElementOfRow + minMaxBoxWidthPoints; columnIndexInRow++)
                {
                    if (terrain.Coordinates[columnIndexInRow].Y >= startPoint.Y && terrain.Coordinates[columnIndexInRow].Y <= endPoint.Y)
                    {
                        if (terrain.Coordinates[columnIndexInRow].Z <= waterLevelHeight)
                        {
                            terrain.Coordinates[columnIndexInRow] = new Vector3(terrain.Coordinates[columnIndexInRow].X, terrain.Coordinates[columnIndexInRow].Y, waterLevelHeight);
                            int indexWidth = columnIndexInRow % NumberOfPointsWidth;
                            int indexHeight = (int)Math.Floor((double)columnIndexInRow / NumberOfPointsWidth);
                            terrain.HeightmapAsArray[indexHeight, indexWidth] = waterLevelHeight;
                           
                        }
                    }
                }
            }
        }

        private List<Edge> CreateTerrainEdgesFromWaterLevelPoints(Terrain terrain)
        {
            List<Edge> edges = new List<Edge>();
            Vector3 start = new Vector3(0, 0, 0);
            if (waterLevelPoints.Count() > 0) start = terrain.Coordinates[GetWaterLevelCoordinateIndex(waterLevelPoints[0])];
            for (int i = 1; i < waterLevelPoints.Count(); i++)
            {
                Vector3 end = terrain.Coordinates[GetWaterLevelCoordinateIndex(waterLevelPoints[i])];
                if (start.X == end.X)
                {
                    if (end.Y >= start.Y)
                    {
                        edges.Add(new Edge(start, end));
                    }
                    else
                    {
                        edges.Add(new Edge(end, start));
                    }
                }
                else if (start.X < end.X)
                {
                    edges.Add(new Edge(start, end));  // vsetky usecky ukladame aby prva suradnica bola mensia
                }
                else
                {
                    edges.Add(new Edge(end, start));
                }
                start = end;
            }
            if (waterLevelPoints.Count() > 0)
            {
                if (start.X < terrain.Coordinates[GetWaterLevelCoordinateIndex(waterLevelPoints[0])].X)
                {
                    edges.Add(new Edge(start, terrain.Coordinates[GetWaterLevelCoordinateIndex(waterLevelPoints[0])]));
                }
                else if (start.X > terrain.Coordinates[GetWaterLevelCoordinateIndex(waterLevelPoints[0])].X)
                {
                    edges.Add(new Edge(terrain.Coordinates[GetWaterLevelCoordinateIndex(waterLevelPoints[0])], start));
                }
                else
                {
                    if (terrain.Coordinates[GetWaterLevelCoordinateIndex(waterLevelPoints[0])].Y >= start.Y)
                    {
                        edges.Add(new Edge(start, terrain.Coordinates[GetWaterLevelCoordinateIndex(waterLevelPoints[0])]));
                    }
                    else
                    {
                        edges.Add(new Edge(terrain.Coordinates[GetWaterLevelCoordinateIndex(waterLevelPoints[0])], start));
                    }
                }
            }
            return edges;
        }

        private void ScanLine(Terrain terrain)
        {
            List<Edge> edges = CreateTerrainEdgesFromWaterLevelPoints(terrain);

            Dictionary<double, List<PolyEdge>> TH = new Dictionary<double, List<PolyEdge>>();
            List<PolyEdge> TAH = new List<PolyEdge>();
            List<PolyEdge> TAHhelpList = new List<PolyEdge>();
            List<Edge> edgesHelpList = new List<Edge>();

            List<System.Windows.Point> minMaxBox = GetMinMaxBoxOfWaterLevel(waterLevelPoints);
            int startIndexOfBox = GetWaterLevelCoordinateIndex(minMaxBox[0]);
            int endIndexOfBox = GetWaterLevelCoordinateIndex(minMaxBox[1]);
            float distanceBetweenPoints = 1;
            if (startIndexOfBox + NumberOfPointsWidth < terrain.Coordinates.Count()) distanceBetweenPoints = terrain.Coordinates[startIndexOfBox + NumberOfPointsWidth].X - terrain.Coordinates[startIndexOfBox].X;

            List<int> indicesOfParallelEdges = new List<int>();
            for (int u = 0; u < edges.Count(); u++)
            {
                if (edges[u].start.X == edges[u].end.X)
                {
                    indicesOfParallelEdges.Add(u);
                }
            }

            // shortening of extremal points
            for (int u = 0; u < edges.Count(); u++)
            {
                Edge nextEdge = edges[0];
                foreach (Edge e in edges)
                {
                    if (edges[u] != e)
                    {                        
                        if (edges[u].end == e.start || edges[u].end == e.end)
                        {
                            nextEdge = e;
                        }
                    }
                }

                Vector3 nextPointToEdge = nextEdge.end;
                if (edges[u].end == nextEdge.end) nextPointToEdge = nextEdge.start;
                if (nextPointToEdge.X >= edges[u].end.X && edges[u].end.X != edges[u].start.X)
                {
                    edges[u].end.X = terrain.Coordinates[terrain.Coordinates.IndexOf(edges[u].end) - NumberOfPointsWidth].X;
                    float smer = (edges[u].end.Y - edges[u].start.Y) / (edges[u].end.X - edges[u].start.X);
                    edges[u].end.Y = edges[u].end.Y - smer * distanceBetweenPoints;
                }
            }


            // remove the horizontal edges
            edgesHelpList = edges.ToList();    
            foreach (int i in indicesOfParallelEdges)
            {
                SetHeigntOnInterval(new Vector3(edgesHelpList[i].start.X, edgesHelpList[i].start.Y, 0), new Vector3(edgesHelpList[i].start.X, edgesHelpList[i].end.Y, 0), terrain);
                edges.Remove(edgesHelpList[i]);
            }

            for (int firstElementOfRow = startIndexOfBox; firstElementOfRow <= endIndexOfBox; firstElementOfRow += NumberOfPointsWidth)
            {
                float i = terrain.Coordinates[firstElementOfRow].X;
                TH.Add(i, new List<PolyEdge>());  // initialize the table of edges              

                foreach (Edge e in edges)
                {
                    if (e.start.X == i)
                    {
                        float directOfEdge = (e.end.Y - e.start.Y) / (e.end.X - e.start.X);
                        TH[i].Add(new PolyEdge(e.end.X, e.start.Y, directOfEdge));
                    }
                }

            }

            for (int firstElementOfRow = startIndexOfBox; firstElementOfRow <= endIndexOfBox; firstElementOfRow += NumberOfPointsWidth)
            {
                float i = terrain.Coordinates[firstElementOfRow].X;
                foreach (PolyEdge k in TH[i])
                {
                    TAH.Add(k);  // active edges on currelt line           
                }

                TAH = TAH.OrderBy(s => s.y).ToList();
                for (int k = 0; k < TAH.Count(); k = k + 2)
                {
                    SetHeigntOnInterval(new Vector3(i, TAH[k].y, 0), new Vector3(i, TAH[k + 1].y, 0), terrain);
                }

                TAHhelpList = TAH.ToList();  
                foreach (PolyEdge activeEdges in TAH)
                {
                    if (activeEdges.xmax == i)
                    {
                        TAHhelpList.Remove(activeEdges);    // remove inactive edges
                    }
                    activeEdges.y = activeEdges.y + activeEdges.direct * distanceBetweenPoints;
                }
                TAH = TAHhelpList.ToList();
            }

        }

        public class PolyEdge
        {   
            public float y;
            public float xmax;
            public float direct;

            public PolyEdge(float xmax, float y, float direct)
            {
                this.xmax = xmax;
                this.y = y;
                this.direct = direct;
            }
        }

        public class Edge
        {  
            public Vector3 start;
            public Vector3 end;

            public Edge(Vector3 start, Vector3 end)
            {
                this.start = start;
                this.end = end;
            }
        }
                
    }
}
