using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TVMS
{
    class Program
    {
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

        private static void Main(string[] args)
        {
            int n = 3900, m = 10;
            var matrix = ReadDataMatrix(@"D:\Универ\3 курс\2 семестр\Тер вер и мат стат\movie_metadata_1.csv");
            var transpMatrix = MatrixFunction.TransposeMatrix(matrix);
            Dictionary<double, double> laplasMatrix = ReadLaplasMatrix("2.txt");

            for (int k = 0; k < m; ++k)
            {
                PearsonConsentCriterion pcc = new PearsonConsentCriterion(transpMatrix[k], laplasMatrix);
                KolmogorovConsentCriterion kcc = new KolmogorovConsentCriterion(transpMatrix[k]);

                Console.WriteLine("{0,20}  {1,20}        {2,15}{3,15}", "     ", "     ", "     Частота", " Вероятность");
                for (int j = 0; j < pcc.Intervals.Length - 1; ++j)
                    Console.WriteLine("{0,20}   -   {1,20} {2,15}{3,15}", 
                        pcc.Intervals[j], pcc.Intervals[j + 1], pcc.HitsInIntervalsCount[j], pcc.HitsInIntervalsProbability[j]);

                Console.WriteLine();
                Console.WriteLine("Среднее значение x: " + pcc.AverageValueX);
                Console.WriteLine("Среднее квадратичное отклонение: " + pcc.MeanSquareDeviation);
                Console.WriteLine("Критерий Пирсона: " + pcc.PearsonCriterionValue);
                Console.WriteLine("Табличный критерий Пирсона = " + 20.519);
                Console.WriteLine("Критерий Колмогорова = " + kcc.KolmogorovCriterionValue);
                Console.WriteLine("Табличный критерий Колмогорова = 1,950");
                Console.WriteLine();
                Console.WriteLine();
            }

            CorrelationsAnalysis correlationsAnalyses = new CorrelationsAnalysis(transpMatrix, m - 1);
            Console.WriteLine("Матрица корреляции");
            for (int u = 0; u < correlationsAnalyses.PairCorrelationsMatrix.Length; u++)
            {
                for (int j = 0; j < correlationsAnalyses.PairCorrelationsMatrix[u].Length; j++)
                    Console.Write("{0}         ", Math.Round(correlationsAnalyses.PairCorrelationsMatrix[u][j], 5));
                Console.WriteLine();
            }

            Console.WriteLine(); Console.WriteLine("Коэффициенты значимости для матрицы парных корреляций: ");
            for (int u = 0; u < correlationsAnalyses.MatrixSignificanceFactors.Length; u++)
            {
                for (int j = 0; j < correlationsAnalyses.MatrixSignificanceFactors[u].Length; j++)
                    Console.Write("{0}         ", Math.Round(correlationsAnalyses.MatrixSignificanceFactors[u][j], 5));
                Console.WriteLine();
            }
            Console.WriteLine();

            Console.WriteLine("Коэффициенты значимости:");
            for (int u = 0; u < correlationsAnalyses.ParametersSignificanceFactors.Length; u++)
                Console.WriteLine(correlationsAnalyses.ParametersSignificanceFactors[u]);
            Console.WriteLine();

            var multipleCoefficientY = correlationsAnalyses.SelectiveMultipleCoefficient;
            var determinationCoefficientY = correlationsAnalyses.DeterminationCoefficient;
            Console.WriteLine("Выборочный множественный коэффициент Y: " + multipleCoefficientY);
            Console.WriteLine("Коэффициент детерминации: " + determinationCoefficientY);
            if (determinationCoefficientY > 0.8)
                Console.WriteLine("Модель адекватна");

            Console.WriteLine();
            Console.WriteLine("Матрица частной корреляции");
            for (int u = 0; u < correlationsAnalyses.PartialCorrelationsMatrix.Length; u++)
            {
                for (int j = 0; j < correlationsAnalyses.PartialCorrelationsMatrix.Length; j++)
                    Console.Write("{0,6}  ", Math.Round(correlationsAnalyses.PartialCorrelationsMatrix[u][j], 4));
                Console.WriteLine();
            }
            Console.WriteLine();

            RegressionAnalysis  ra = new RegressionAnalysis(matrix);
            //Console.WriteLine("Коэффициенты регрессии:");
            //for (int i = 0; i < ra.RegressionCoefficients.Length; i++)
            //{
            //    if (i == 0)
            //        Console.WriteLine("a = {0}", ra.RegressionCoefficients[i]);
            //    else
            //        Console.WriteLine("b{0} = {1}", i, ra.RegressionCoefficients[i]);
            //}

            Console.WriteLine("Значимость коэффициентов регрессии:");
            for (int i = 0; i < m; i++)
                Console.WriteLine("b{0} = {1}", (i + 1), ra.RegressionCoefficientsSignificance[i]);

            Console.WriteLine("Доверительные интервалы коэффициентов регрессии:");
            for (int i = 0; i < m; i++)
                Console.WriteLine("{0} <= b{1} <= {2}", ra.ConfidenceIntervalsOfCoefficients[i].Item1, (i + 1), ra.ConfidenceIntervalsOfCoefficients[i].Item2);

            Console.WriteLine("Коэффициент значимости уравнения регрессии:");
            Console.WriteLine(ra.RegressionEquationSignificance);
            Console.WriteLine();
            Console.WriteLine("Коэффициенты эластичности");
            double[] elast = ra.ElasticityCoefficients;
            for (int u = 0; u < m; u++)
                Console.WriteLine("x{0} = {1}", (u + 1), elast[u]);

            //double[] X0 = { 1, 0.92, 0.83, 5.82, 1.2, 4.25, 1.01, 1.01, 1 };
            //double[] resultat = Prognoz(matrix, X0);

            Forecast f = new Forecast(matrix);
            Console.WriteLine("Прогнозирование");
            Console.WriteLine(f.Value[0] + " < y < " + f.Value[1]);
            Console.WriteLine();
        }
    }
}