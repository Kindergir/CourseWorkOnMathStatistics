using System;
using System.Linq;

namespace MathStatisticsLibrary
{
    public class KolmogorovConsentCriterion
    {
        public double KolmogorovCriterionValue { get; }

        public KolmogorovConsentCriterion(double[] parameterValues)
        {
            KolmogorovCriterionValue = (6 * Math.Sqrt(parameterValues.Length) * KolmogorovTest(parameterValues) + 1) 
                / (6 * Math.Sqrt(parameterValues.Length));
        }

        private double KolmogorovTest(double[] parameterValues)
        {
            int n = parameterValues.Length;
            double[] d1 = new double[parameterValues.Length];
            double[] d2 = new double[parameterValues.Length];

            for (int i = 0; i < parameterValues.Length; ++i)
                d1[i] = ((double)i / n) - (Math.Pow(Math.E, (-Math.Pow(parameterValues[i], 2) / 2) / Math.Sqrt(2 * Math.PI)));

            for (int i = 0; i < parameterValues.Length; ++i)
                d2[i] = (Math.Pow(Math.E, (-Math.Pow(parameterValues[i], 2) / 2) / Math.Sqrt(2 * Math.PI))) - ((double)i - 1) / n;

            return Math.Max(d1.Max(), d2.Max());
        }
    }
}
