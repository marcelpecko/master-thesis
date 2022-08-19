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
using System.Windows.Shapes;

namespace ProceduralTerrainGeneration
{   
    public partial class StartWindow : Window
    {
        public StartWindow()
        {
            InitializeComponent();
        }

        private void SelectButton(object sender, RoutedEventArgs e)
        {
            string terrainType = (sender as Button).Name;
            MainWindow mainWindow = new MainWindow(terrainType);
            mainWindow.Show();
            Close();
        }
    }
}
