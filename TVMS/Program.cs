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
        static double SelectiveMultipleCoefficient(double[][] r, int i) 
        {
            return Math.Sqrt(1 - (MatrixFunction.MatrixDeterminant(r) / MatrixFunction.AlgebraicComplement(r, i, i)));
        }

        public static double[] Regr(double[][] A)
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

            result = MatrixFunction.MatrixVectorMultiplication(MatrixFunction.MultiplicateMatrix(MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(A1), A1)), MatrixFunction.TransposeMatrix(A1)), Y);

            double Qr = MatrixFunction.SumOfVectorCoordinatesSquares(MatrixFunction.MatrixVectorMultiplication(A1, result));
            double Qost = MatrixFunction.SumOfVectorCoordinatesSquares(MatrixFunction.VectorsDifference(Y, MatrixFunction.MatrixVectorMultiplication(A1, result)));

            double s = Qost / (n - m - 1);
            double[] sb = new double[result.Length];
            var temperory = (MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(A1), A1)));

            for (int j = 0; j < m; j++)
                sb[j] = s * temperory[j][j];

            double t = 4.587;
            Console.WriteLine("Коэффициенты регрессии:");

            for (int i = 0; i < result.Length; i++)
            {
                if (i == 0)
                    Console.WriteLine("a = {0}", result[i]);
                else
                    Console.WriteLine("b{0} = {1}", i, result[i]);
            }

            Console.WriteLine("Значимость коэффициентов регрессии:");
            for (int i = 0; i < m; i++)
                Console.WriteLine("b{0} = {1}", (i + 1), result[i] / sb[i]);

            Console.WriteLine("Доверительные интервалы коэффициентов регрессии:");
            for (int i = 0; i < m; i++)
                Console.WriteLine("{0} <= b{1} <= {2}", (result[i] - t * sb[i]), (i + 1), (result[i] + t * sb[i]));

            Console.WriteLine("Коэффициент значимости уравнения регрессии:");
            Console.WriteLine(Qr * (n - m - 1) / (Qost * (m + 1)));
            return result;
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

            PairCorrelationsMatrix pcm = new PairCorrelationsMatrix(transpMatrix);
            Console.WriteLine("Матрица корреляции");
            for (int u = 0; u < pcm.CorrelationsMatrix.Length; u++)
            {
                for (int j = 0; j < pcm.CorrelationsMatrix[u].Length; j++)
                    Console.Write("{0}         ", Math.Round(pcm.CorrelationsMatrix[u][j], 5));
                Console.WriteLine();
            }

            Console.WriteLine(); Console.WriteLine("Коэффициенты значимости для матрицы парных корреляций: ");
            for (int u = 0; u < pcm.SignificanceFactorsMatrix.Length; u++)
            {
                for (int j = 0; j < pcm.SignificanceFactorsMatrix[u].Length; j++)
                    Console.Write("{0}         ", Math.Round(pcm.SignificanceFactorsMatrix[u][j], 5));
                Console.WriteLine();
            }
            Console.WriteLine();

            Console.WriteLine("Коэффициенты значимости:");
            for (int u = 0; u < pcm.CorrelationsMatrix.Length; u++)
                Console.WriteLine((Math.Pow(SelectiveMultipleCoefficient(pcm.CorrelationsMatrix, u), 2)) / ((1 - Math.Pow(SelectiveMultipleCoefficient(pcm.CorrelationsMatrix, u), 2)) / (n - 2)));
            Console.WriteLine();

            double Rmn1 = SelectiveMultipleCoefficient(pcm.CorrelationsMatrix, 1);
            double D = Math.Pow(SelectiveMultipleCoefficient(pcm.CorrelationsMatrix, 1), 2);
            Console.WriteLine("Выборочный множественный коэффициент Y: " + Rmn1);
            Console.WriteLine("Коэффициент детерминации: " + D);
            if (D > 0.8)
                Console.WriteLine("Модель адекватна");

            Console.WriteLine();
            Console.WriteLine("Матрица частной корреляции");
            double[][] r1 = PartialCorrelationsMatrix.GetForMatrix(pcm.CorrelationsMatrix);
            for (int u = 0; u < r1.Length; u++)
            {
                for (int j = 0; j < r1.Length; j++)
                    Console.Write("{0,6}  ", Math.Round(r1[u][j], 4));
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