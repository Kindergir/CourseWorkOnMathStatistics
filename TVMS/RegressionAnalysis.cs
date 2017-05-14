using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TVMS
{
    public class RegressionAnalysis
    {
        public double[] RegressionCoefficients { get; }
        public double[] RegressionCoefficientsSignificance { get; }
        public double RegressionEquationSignificance { get; }
        public Tuple<double, double> ConfidenceIntervalsOfCoefficients { get; }
    }
}
