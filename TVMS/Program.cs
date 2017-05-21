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

        private static void WriteMatrix(string path, double[][] matrix)
        {
            StreamWriter sw = new StreamWriter(path);
            for (int i = 0; i < matrix.Length; ++i)
            {
                string cur = "";
                for (int j = 0; j < matrix[i].Length; ++j)
                {
                    cur += matrix[i][j] + ";";
                }
                if (cur[cur.Length - 1] == ';')
                    cur = cur.Substring(0, cur.Length - 1);
                cur += "\n";
                sw.Write(cur);
            }
        }

        private static void Main(string[] args)
        {
            int n = 629, m = 10;
            var matrix = ReadDataMatrix(@"D:\Универ\3 курс\2 семестр\Тер вер и мат стат\movie_metadata_1.csv");
            var transpMatrix = MatrixFunction.TransposeMatrix(matrix);
            Dictionary<double, double> laplasMatrix = ReadLaplasMatrix("2.txt");

            for (int k = 0; k < m; ++k)
            {
                PearsonConsentCriterion pcc = new PearsonConsentCriterion(transpMatrix[k].Take(n).ToArray(), laplasMatrix);
                KolmogorovConsentCriterion kcc = new KolmogorovConsentCriterion(transpMatrix[k]);
                DescriptiveStatistics ds = new DescriptiveStatistics(transpMatrix[k]);

                Console.WriteLine();
                Console.WriteLine("Среднее арифметическое = {0}", ds.ArithmeticalMean);
                Console.WriteLine("Мода = {0}", ds.Mode);
                Console.WriteLine("Медиана = {0}", ds.Median);
                Console.WriteLine("Дисперсия = {0}", ds.Dispersion);
                Console.WriteLine("Асимметрия = {0}", ds.Assimmetry);
                Console.WriteLine("Эксцесс = {0}", ds.Excess);
                Console.WriteLine("Стандартное отклонение = {0}", ds.StandardDeviation);
                Console.WriteLine("Коэффициент вариации = {0}", ds.VariationCoefficient);
                Console.WriteLine("Размах вариации = {0}", ds.VariationRange);
                Console.WriteLine("Среднее значение x: " + pcc.AverageValueX);
                Console.WriteLine("Среднее квадратичное отклонение: " + pcc.MeanSquareDeviation);
                Console.WriteLine("Критерий Пирсона: " + pcc.PearsonCriterionValue);
                Console.WriteLine("Табличный критерий Пирсона = " + 63.6567);
                Console.WriteLine("Критерий Колмогорова = " + kcc.KolmogorovCriterionValue);
                Console.WriteLine("Табличный критерий Колмогорова = 1,950\n\n");
                Console.WriteLine();
            }

            CorrelationsAnalysis correlationsAnalyses = new CorrelationsAnalysis(transpMatrix, m - 1);
            //WriteMatrix(@"D:\Универ\3 курс\2 семестр\Тер вер и мат стат\PairCorrelationMatrix.csv", correlationsAnalyses.PairCorrelationsMatrix);
            Console.WriteLine("Матрица корреляции");
            foreach (double[] correlationCoef in correlationsAnalyses.PairCorrelationsMatrix)
            {
                for (int j = 0; j < correlationCoef.Length; j++)
                    Console.Write("{0}    ", Math.Round(correlationCoef[j], 5));
                Console.WriteLine();
            }

            Console.WriteLine(); Console.WriteLine("Коэффициенты значимости для матрицы парных корреляций: ");
            foreach (double[] t in correlationsAnalyses.MatrixSignificanceFactors)
            {
                for (int j = 0; j < t.Length; j++)
                    Console.Write("{0}         ", Math.Round(t[j], 5));
                Console.WriteLine();
            }
            Console.WriteLine();

            Console.WriteLine("Коэффициенты значимости:");
            foreach (double t in correlationsAnalyses.ParametersSignificanceFactors)
                Console.WriteLine(t);
            Console.WriteLine();

            var multipleCoefficientY = correlationsAnalyses.SelectiveMultipleCoefficient;
            var determinationCoefficientY = correlationsAnalyses.DeterminationCoefficient;
            Console.WriteLine("Выборочный множественный коэффициент Y: " + multipleCoefficientY);
            Console.WriteLine("Коэффициент детерминации: " + determinationCoefficientY);
            if (determinationCoefficientY > 0.8)
                Console.WriteLine("Модель адекватна");

            Console.WriteLine();
            Console.WriteLine("Матрица частной корреляции");
            foreach (double[] t in correlationsAnalyses.PartialCorrelationsMatrix)
            {
                for (int j = 0; j < correlationsAnalyses.PartialCorrelationsMatrix.Length; j++)
                    Console.Write("{0,6}  ", Math.Round(t[j], 4));
                Console.WriteLine();
            }
            Console.WriteLine();

            RegressionAnalysis ra = new RegressionAnalysis(matrix, m - 1);
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