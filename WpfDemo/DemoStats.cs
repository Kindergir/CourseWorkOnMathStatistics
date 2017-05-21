using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathStatisticsLibrary;

namespace WpfDemo
{
    public class DemoStats
    {
        public CorrelationsAnalysis Correlations { get; private set; }
        public RegressionAnalysis Regressions { get; private set; }
        public Forecast Forecast { get; private set; }
        public DescriptiveStatistics[] Descriptive { get; private set; }

        public DemoStats(double[][] dataMatrix)
        {
            var transpMatrix = MatrixFunction.TransposeMatrix(dataMatrix);
            Correlations = new CorrelationsAnalysis(transpMatrix, dataMatrix[0].Length - 1);
            Regressions = new RegressionAnalysis(dataMatrix, dataMatrix[0].Length - 1);
            Descriptive = new DescriptiveStatistics[transpMatrix.Length];
            for (int i = 0; i < Descriptive.Length; ++i)
                Descriptive[i] = new DescriptiveStatistics(transpMatrix[i]);
            Forecast = new Forecast(dataMatrix, transpMatrix.Length - 1);
        }
    }
}
