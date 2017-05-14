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
        public static double[] Regr(double[][] dataMatrix)
        {
            int n = dataMatrix.Length, parametersCount = dataMatrix[0].Length;
            double[] regressionCoefficients = new double[parametersCount];
            double[] resultParameter = new double[n];
            double[][] A1 = new double[n][];

            for (int i = 0; i < n; i++)
            {
                A1[i] = new double[parametersCount];
                resultParameter[i] = dataMatrix[i][parametersCount - 1];
            }

            for (int i = 0; i < n; i++)
                A1[i][parametersCount - 1] = 1;

            for (int i = 0; i < n; i++)
                for (int j = 0; j < parametersCount - 1; j++)
                    A1[i][j] = dataMatrix[i][j];

            regressionCoefficients = MatrixFunction.MatrixVectorMultiplication(MatrixFunction.MultiplicateMatrix(MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(A1), A1)), MatrixFunction.TransposeMatrix(A1)), resultParameter);

            double Qr = MatrixFunction.SumOfVectorCoordinatesSquares(MatrixFunction.MatrixVectorMultiplication(MatrixFunction.TransposeMatrix(A1), regressionCoefficients));
            double Qost = MatrixFunction.SumOfVectorCoordinatesSquares(MatrixFunction.VectorsDifference(resultParameter, MatrixFunction.MatrixVectorMultiplication(MatrixFunction.TransposeMatrix(A1), regressionCoefficients)));

            double s = Qost / (n - parametersCount - 1);
            double[] sb = new double[parametersCount];
            var temperory = MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(A1), A1));

            for (int j = 0; j < parametersCount; j++)
                sb[j] = s * temperory[j][j];

            double t = 4.587;
            Console.WriteLine("Коэффициенты регрессии:");
            for (int i = 0; i < regressionCoefficients.Length; i++)
            {
                if (i == 0)
                    Console.WriteLine("a = {0}", regressionCoefficients[i]);
                else
                    Console.WriteLine("b{0} = {1}", i, regressionCoefficients[i]);
            }

            Console.WriteLine("Значимость коэффициентов регрессии:");
            for (int i = 0; i < parametersCount; i++)
                Console.WriteLine("b{0} = {1}", (i + 1), regressionCoefficients[i] / sb[i]);

            Console.WriteLine("Доверительные интервалы коэффициентов регрессии:");
            for (int i = 0; i < parametersCount; i++)
                Console.WriteLine("{0} <= b{1} <= {2}", (regressionCoefficients[i] - t * sb[i]), (i + 1), (regressionCoefficients[i] + t * sb[i]));

            Console.WriteLine("Коэффициент значимости уравнения регрессии:");
            Console.WriteLine(Qr * (n - parametersCount - 1) / (Qost * (parametersCount + 1)));
            return regressionCoefficients;
        }

        public static double[] ElasticityCoefficient(double[][] A, double[] b)
        {
            double[] result = new double[A.Length];

            double y = A[0].Average();
            result[0] = b[0];
            for (int i = 1; i < A.Length; ++i)
                result[i] = b[i] * A[i].Average() / y;

            return result;
        }

        public static double[] Prognoz(double[][] A, double[] X0)
        {
            int n = A.Length, m = A[0].Length;
            double[] result = new double[m];
            double[] Y = new double[n];
            double[][] A1 = new double[n][];
            for (int i = 0; i < n; i++)
            {
                A1[i] = new double[m];
                Y[i] = A[i][0];
            }

            for (int i = 0; i < n; i++)
                A1[i][0] = 1;

            for (int i = 0; i < n; i++)
                for (int j = 1; j < m; j++)
                    A1[i][j] = A[i][j];

            result = MatrixFunction.MatrixVectorMultiplication((MatrixFunction.MultiplicateMatrix(MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(A1), A1)), MatrixFunction.TransposeMatrix(A1))), Y);
            double Qost = MatrixFunction.SumOfVectorCoordinatesSquares(MatrixFunction.VectorsDifference(Y, MatrixFunction.MatrixVectorMultiplication(A1, result)));
            double s = Math.Sqrt(Qost / (n - m - 1));
            double t = 4.587;
            double[] a = new double[2];
            a[0] = MatrixFunction.ScalarProductOfVectors(X0, result) - t * s * Math.Sqrt(MatrixFunction.ScalarProductOfVectors(MatrixFunction.TransposeMatrixVectorProduct(X0, MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(A1), A1))), X0));
            a[1] = MatrixFunction.ScalarProductOfVectors(X0, result) + t * s * Math.Sqrt(MatrixFunction.ScalarProductOfVectors(MatrixFunction.TransposeMatrixVectorProduct(X0, MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(A1), A1))), X0));
            return a;
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

            double[] result = Regr(matrix);
            Console.WriteLine();
            Console.WriteLine("Коэффициенты эластичности");
            double[] elast = ElasticityCoefficient(transpMatrix, result);
            for (int u = 0; u < m; u++)
                Console.WriteLine("x{0} = {1}", (u + 1), elast[u]);

            double[] X0 = { 1, 0.92, 0.83, 5.82, 1.2, 4.25, 1.01, 1.01, 1 };
            double[] resultat = Prognoz(matrix, X0);
            Console.WriteLine("Прогнозирование");
            Console.WriteLine(resultat[0] + " < y < " + resultat[1]);
            Console.WriteLine();
        }
    }
}