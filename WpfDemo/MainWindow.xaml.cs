using System;
using System.Collections.Generic;
using System.IO;
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
using MathStatisticsLibrary;
using Microsoft.Win32;

namespace WpfDemo
{
    /// <summary>
    /// Логика взаимодействия для MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();

        }

        private void BtnCalc_OnClick(object sender, RoutedEventArgs e)
        {
            int n = 629, m = 10;
            var matrix = ReadDataMatrix(@"D:\Универ\3 курс\2 семестр\Тер вер и мат стат\movie_metadata_1.csv");
            var transpMatrix = MatrixFunction.TransposeMatrix(matrix);
            Dictionary<double, double> laplasMatrix = ReadLaplasMatrix("2.txt");

            tbDescript.Text = "";
            for (int k = 0; k < m; ++k)
            {
                PearsonConsentCriterion pcc = new PearsonConsentCriterion(transpMatrix[k].Take(n).ToArray(), laplasMatrix);
                KolmogorovConsentCriterion kcc = new KolmogorovConsentCriterion(transpMatrix[k]);
                DescriptiveStatistics ds = new DescriptiveStatistics(transpMatrix[k]);
                tbDescript.Text += String.Format("Параметр: {0}\n", k + 1);
                tbDescript.Text += String.Format("Среднее арифметическое = {0}\n", ds.ArithmeticalMean);
                tbDescript.Text += String.Format("Мода = {0}\n", ds.Mode);
                tbDescript.Text += String.Format("Медиана = {0}\n", ds.Median);
                tbDescript.Text += String.Format("Дисперсия = {0}\n", ds.Dispersion);
                tbDescript.Text += String.Format("Асимметрия = {0}\n", ds.Assimmetry);
                tbDescript.Text += String.Format("Эксцесс = {0}\n", ds.Excess);
                tbDescript.Text += String.Format("Стандартное отклонение = {0}\n", ds.StandardDeviation);
                tbDescript.Text += String.Format("Коэффициент вариации = {0}\n", ds.VariationCoefficient);
                tbDescript.Text += String.Format("Размах вариации = {0}\n", ds.VariationRange);
                tbDescript.Text += String.Format("Среднее значение x: " + pcc.AverageValueX);
                tbDescript.Text += String.Format("Среднее квадратичное отклонение: " + pcc.MeanSquareDeviation);
                tbDescript.Text += String.Format("Критерий Пирсона: " + pcc.PearsonCriterionValue);
                tbDescript.Text += String.Format("Табличный критерий Пирсона = " + 63.6567);
                tbDescript.Text += String.Format("Критерий Колмогорова = " + kcc.KolmogorovCriterionValue);
                tbDescript.Text += String.Format("Табличный критерий Колмогорова = 1,950\n\n");
            }

            tbCorrelations.Text = "";
            CorrelationsAnalysis correlationsAnalyses = new CorrelationsAnalysis(transpMatrix, m - 1);
            tbCorrelations.Text += "Матрица корреляции\n";
            foreach (double[] correlationCoef in correlationsAnalyses.PairCorrelationsMatrix)
            {
                for (int j = 0; j < correlationCoef.Length; j++)
                    tbCorrelations.Text += String.Format("{0}    ", Math.Round(correlationCoef[j], 5));
                tbCorrelations.Text += '\n';
            }

            tbCorrelations.Text += "Коэффициенты значимости для матрицы парных корреляций: \n";
            foreach (double[] t in correlationsAnalyses.MatrixSignificanceFactors)
            {
                for (int j = 0; j < t.Length; j++)
                    tbCorrelations.Text += String.Format("{0}         ", Math.Round(t[j], 5));
                tbCorrelations.Text += '\n';
            }

            tbCorrelations.Text += "Коэффициенты значимости:\n";
            foreach (double t in correlationsAnalyses.ParametersSignificanceFactors)
                tbCorrelations.Text += t;
            tbCorrelations.Text += "\n";

            var multipleCoefficientY = correlationsAnalyses.SelectiveMultipleCoefficient;
            var determinationCoefficientY = correlationsAnalyses.DeterminationCoefficient;
            tbCorrelations.Text += String.Format("Выборочный множественный коэффициент Y: " + multipleCoefficientY);
            tbCorrelations.Text += String.Format("Коэффициент детерминации: " + determinationCoefficientY);
            if (determinationCoefficientY > 0.75)
                tbCorrelations.Text += String.Format("Модель адекватна");

            tbCorrelations.Text += '\n';
            tbCorrelations.Text += String.Format("Матрица частной корреляции");
            foreach (double[] t in correlationsAnalyses.PartialCorrelationsMatrix)
            {
                for (int j = 0; j < correlationsAnalyses.PartialCorrelationsMatrix.Length; j++)
                    tbCorrelations.Text += String.Format("{0,6}  ", Math.Round(t[j], 4));
                tbCorrelations.Text += "\n";
            }
            tbCorrelations.Text += "\n";

            RegressionAnalysis ra = new RegressionAnalysis(matrix, m - 1);
            tbRegressions.Text = "Коэффициенты регрессии:\n";
            for (int i = 0; i < ra.RegressionCoefficients.Length; i++)
            {
                if (i == 0)
                    tbRegressions.Text += String.Format("a = {0}\n", ra.RegressionCoefficients[i]);
                else
                    tbRegressions.Text += String.Format("b{0} = {1}\n", i, ra.RegressionCoefficients[i]);
            }

            tbRegressions.Text += String.Format("Значимость коэффициентов регрессии:");
            for (int i = 0; i < m; i++)
                tbRegressions.Text += String.Format("b{0} = {1}", (i + 1), ra.RegressionCoefficientsSignificance[i]);

            tbRegressions.Text += String.Format("Доверительные интервалы коэффициентов регрессии:");
            for (int i = 0; i < m; i++)
                tbRegressions.Text += String.Format("{0} <= b{1} <= {2}", ra.ConfidenceIntervalsOfCoefficients[i].Item1, (i + 1), ra.ConfidenceIntervalsOfCoefficients[i].Item2);

            tbRegressions.Text += String.Format("Коэффициент значимости уравнения регрессии:");
            tbRegressions.Text += String.Format(ra.RegressionEquationSignificance.ToString());

            tbRegressions.Text += String.Format("\nКоэффициенты эластичности");
            double[] elast = ra.ElasticityCoefficients;
            for (int u = 0; u < m; u++)
                tbRegressions.Text += String.Format("x{0} = {1}", (u + 1), elast[u]);

            Forecast f = new Forecast(matrix, m - 1);
            tbForecast.Text = "Прогнозирование\n";
            tbForecast.Text += string.Format(f.Value[0] + " < y < " + f.Value[1]);
        }

        private void WriteCorrelationsMatrix(double[][] matrix)
        {
            tbCorrelations.Text = "";
            foreach (double[] correlationCoef in matrix)
            {
                for (int j = 0; j < correlationCoef.Length; j++)
                    tbCorrelations.Text += String.Format("{0}    ", Math.Round(correlationCoef[j], 5));
                tbCorrelations.Text += '\n';
            }
        }

        private static double[][] ReadDataMatrix(string path)
        {
            int n = 3900, m = 10;

            double[][] matrix = new double[n][];
            for (int k = 0; k < n; k++)
                matrix[k] = new double[m];

            StreamReader fn = new StreamReader(path, Encoding.Default);
            for (int i = 0; !fn.EndOfStream; i++)
            {
                string[] s = fn.ReadLine().Split(';');
                for (int j = 0; j < m; j++)
                    try
                    {
                        matrix[i][j] = double.Parse(s[j]);
                    }
                    catch { }
            }
            fn.Close();
            return matrix;
        }

        private static Dictionary<double, double> ReadLaplasMatrix(string path)
        {
            StreamReader fs = new StreamReader(path, Encoding.Default);
            Dictionary<double, double> laplas = new Dictionary<double, double>();
            while (!fs.EndOfStream)
            {
                string l = fs.ReadLine();
                string[] s = l.Split('\t');
                double[] s1 = new double[s.Length];
                for (int j = 0; j < s.Length; j++)
                    s1[j] = double.Parse(s[j].Replace('.', ','));
                laplas.Add(s1[0], s1[1] / 2);
            }
            fs.Close();
            return laplas;
        }

        private void BtnChooseFile_OnClick(object sender, RoutedEventArgs e)
        {
            OpenFileDialog dialog = new OpenFileDialog();
            dialog.ShowDialog();
            tbMatrixFilePath.Text = dialog.FileName;
        }
    }
}
