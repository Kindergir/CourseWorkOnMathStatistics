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
        public double ExplicatedDispersion { get; private set; }
        public double ResidualDispersion { get; private set; }
        public double[] RegressionCoefficientsSignificance { get; private set; }
        public double RegressionEquationSignificance { get; private set; }
        public Tuple<double, double>[] ConfidenceIntervalsOfCoefficients { get; private set; }
        public double[] ElasticityCoefficients { get; private set; }

        private readonly double[][] dataMatrix;
        private int resultParameterNumber, parametersCount, measuresCount;
        private double[][] xMatrix;
        private double[] yVector;

        public RegressionAnalysis(double[][] dataMatrix, int resultParameterNumber)
        {
            this.dataMatrix = dataMatrix;
            this.resultParameterNumber = resultParameterNumber;
            parametersCount = dataMatrix[0].Length;
            measuresCount = dataMatrix.Length;
            CalcRegressionCoefficients();
            CalcDispersion();
            CalcElasticityCoefficients();
            CalcRegressionCoefficientsSignificance();
            CalcConfidenceIntervalsOfCoefficients();
            CalcRegressionEquationSignificance();
        }

        private void CalcRegressionCoefficients()
        {
            var transpMatrix = MatrixFunction.TransposeMatrix(dataMatrix);
            int parametersCount = dataMatrix[0].Length, measuresCount = dataMatrix.Length;

            double[] resultParameter = transpMatrix[resultParameterNumber];
            List<double>[] withoutResultParameter = new List<double>[measuresCount];

            for (int i = 0; i < measuresCount; i++)
            {
                withoutResultParameter[i] = new List<double>();
                withoutResultParameter[i].Add(1);

                for (int j = 0; j < parametersCount; ++j)
                    if (j != resultParameterNumber)
                        withoutResultParameter[i].Add(dataMatrix[i][j]);
            }

            xMatrix = withoutResultParameter.Select(o => o.ToArray()).ToArray();
            var transpWithoutMatrix = MatrixFunction.TransposeMatrix(xMatrix);

            RegressionCoefficients = MatrixFunction.MatrixVectorMultiplication(MatrixFunction.MultiplicateMatrix(MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(transpWithoutMatrix, xMatrix)), transpWithoutMatrix), resultParameter);
            RegressionCoefficients = RegressionCoefficients.Take(parametersCount).ToArray();
        }

        private void CalcDispersion()
        {
            var yWithWave = MatrixFunction.MatrixVectorMultiplication(xMatrix, RegressionCoefficients);
            var averageY = yWithWave.Average();
            var devivationsYSum = yWithWave.Select(y => Math.Pow(y - averageY, 2)).Sum();

            ExplicatedDispersion = devivationsYSum / parametersCount;
            ResidualDispersion = devivationsYSum / (measuresCount - parametersCount - 1);
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

        private void CalcRegressionCoefficientsSignificance()
        {
            double[] sb = new double[parametersCount];
            var temperory = MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(xMatrix), xMatrix));

            for (int j = 0; j < parametersCount; j++)
                sb[j] = ResidualDispersion * temperory[j][j];

            RegressionCoefficientsSignificance = new double[parametersCount];
            for (int i = 0; i < parametersCount; ++i)
                RegressionCoefficientsSignificance[i] = RegressionCoefficients[i] / sb[i];
        }

        private void CalcRegressionEquationSignificance()
        {
            RegressionEquationSignificance = ExplicatedDispersion * (parametersCount - measuresCount - 1) / (ResidualDispersion * (measuresCount + 1));
        }

        private void CalcConfidenceIntervalsOfCoefficients()
        {
            double[] sb = new double[parametersCount];
            var temperory = MatrixFunction.InverseMatrix(MatrixFunction.MultiplicateMatrix(MatrixFunction.TransposeMatrix(xMatrix), xMatrix));

            for (int j = 0; j < parametersCount; j++)
                sb[j] = ResidualDispersion * temperory[j][j];

            double t = 4.587;
            ConfidenceIntervalsOfCoefficients = new Tuple<double, double>[parametersCount];
            for (int i = 0; i < parametersCount; ++i)
                ConfidenceIntervalsOfCoefficients[i] = new Tuple<double, double>(RegressionCoefficients[i] - t * sb[i], RegressionCoefficients[i] + t * sb[i]);
        }
    }
}
