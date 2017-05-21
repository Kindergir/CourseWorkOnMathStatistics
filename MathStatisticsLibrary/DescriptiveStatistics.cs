using System;
using System.Linq;

namespace MathStatisticsLibrary
{
    public class DescriptiveStatistics
    {
        public double ArithmeticalMean { get; private set; }
        public double Mode { get; private set; }
        public double Median { get; private set; }
        public double VariationRange { get; private set; }
        public double StandardDeviation { get; private set; }
        public double Dispersion { get; private set; }
        public double Assimmetry { get; private set; }
        public double Excess { get; private set; }
        public double VariationCoefficient { get; private set; }

        private double[] variable;

        public DescriptiveStatistics(double[] variable)
        {
            this.variable = variable;
            CalcArithmeticalMean();
            CalcMode();
            CalcMedian();
            CalcVariationRange();
            CalcDispersion();
            CalcStandardDeviation();
            CalcAssimmetry();
            CalcExcess();
            CalcVariationCoefficient();
        }

        private void CalcArithmeticalMean()
        {
            ArithmeticalMean = variable.Sum() / variable.Length;
        }

        private void CalcMode()
        {
            Mode = variable
                .Select(x => new {Value = x, Repeat = variable.Count(y => y == x)})
                .OrderBy(x => x.Repeat)
                .First()
                .Value;
        }

        private void CalcMedian()
        {
            Median = variable
                .OrderBy(x => x)
                .ElementAt(variable.Length/2);
        }

        private void CalcVariationRange()
        {
            VariationRange = variable.Max() - variable.Min();
        }

        private void CalcDispersion()
        {
            Dispersion = variable
                .Select(x => x * x - ArithmeticalMean * ArithmeticalMean)
                .Sum() 
                / (variable.Length - 1);
        }

        private void CalcStandardDeviation()
        {
            StandardDeviation = Math.Sqrt(Dispersion);
        }

        private void CalcAssimmetry()
        {
            Assimmetry = Dispersion / StandardDeviation;
        }

        private void CalcExcess()
        {
            Excess = Assimmetry - 4;
        }

        private void CalcVariationCoefficient()
        {
            VariationCoefficient = StandardDeviation / ArithmeticalMean * 100;
        }
    }
}
