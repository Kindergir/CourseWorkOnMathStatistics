using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TVMS
{
    public class RegressionAnalysis
    {
        public double[] RegressionCoefficients { get; private set; }
        public double[] RegressionCoefficientsSignificance { get; private set; }
        public double RegressionEquationSignificance { get; private set; }
        public Tuple<double, double>[] ConfidenceIntervalsOfCoefficients { get; private set; }
        public double[] ElasticityCoefficients { get; private set; }

        private readonly double[][] dataMatrix;

        public RegressionAnalysis(double[][] dataMatrix)
        {
            this.dataMatrix = dataMatrix;
            Regr();
            CalcElasticityCoefficients();
        }

        private void Regr()
        {
            int n = dataMatrix.Length, parametersCount = dataMatrix[0].Length;
            RegressionCoefficients = new double[parametersCount];
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

            RegressionCoefficients = MatrixFunction.MatrixVectorMultiplication(MatrixFunction.MultiplicateMatrix(MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(A1), A1)), MatrixFunction.TransposeMatrix(A1)), resultParameter);

            double Qr = MatrixFunction.SumOfVectorCoordinatesSquares(MatrixFunction.MatrixVectorMultiplication(MatrixFunction.TransposeMatrix(A1), RegressionCoefficients));
            double Qost = MatrixFunction.SumOfVectorCoordinatesSquares(MatrixFunction.VectorsDifference(resultParameter, MatrixFunction.MatrixVectorMultiplication(MatrixFunction.TransposeMatrix(A1), RegressionCoefficients)));

            double s = Qost / (n - parametersCount - 1);
            double[] sb = new double[parametersCount];
            var temperory = MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(A1), A1));

            for (int j = 0; j < parametersCount; j++)
                sb[j] = s * temperory[j][j];

            RegressionCoefficientsSignificance = new double[parametersCount];
            for (int i = 0; i < parametersCount; ++i)
                RegressionCoefficientsSignificance[i] = RegressionCoefficients[i] / sb[i];

            double t = 4.587;
            ConfidenceIntervalsOfCoefficients = new Tuple<double, double>[parametersCount];
            for (int i = 0; i < parametersCount; ++i)
                ConfidenceIntervalsOfCoefficients[i] = new Tuple<double, double>(RegressionCoefficients[i] - t * sb[i], RegressionCoefficients[i] + t * sb[i]);

            RegressionEquationSignificance = Qr * (n - parametersCount - 1) / (Qost * (parametersCount + 1));
        }
        private void CalcElasticityCoefficients()
        {
            double[] result = new double[dataMatrix.Length];
            var tm = MatrixFunction.TransposeMatrix(dataMatrix);
            double y = tm[tm.Length - 1].Average();
            result[tm.Length - 1] = RegressionCoefficients[tm.Length - 1];
            for (int i = 0; i < tm.Length - 1; ++i)
                result[i] = RegressionCoefficients[i] * tm[i].Average() / y;

            ElasticityCoefficients = result;
        }
    }
}
