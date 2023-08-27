using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace HermiteSpline
{
    public partial class Form1 : System.Windows.Forms.Form
    {
        double a = 2.0, b = 7.0, nu=Double.PositiveInfinity;
        int linknum = 1, funcnum = 9;
        Spline s;

        Thread evaluating;
        
        Dictionary<TrackBar, Chart> tbtochart = new Dictionary<TrackBar, Chart>();

        public Form1()
        {
            InitializeComponent();

            tbtochart.Add(trackBar1, chart1);
            tbtochart.Add(trackBar2, chart2);
            tbtochart.Add(trackBar3, chart3);
            tbtochart.Add(trackBar4, chart4);
            comboBoxLinks.SelectedIndex = linknum - 1;
            comboBoxFunctions.SelectedIndex = funcnum;
        }

        private void textBoxBorderA_Leave(object sender, EventArgs e)
        {
            if (Double.TryParse(((TextBox)sender).Text, out double tmp))
                a = tmp;
            else
                ((TextBox)sender).Text = a.ToString();
        }

        private void textBoxBorderB_Leave(object sender, EventArgs e)
        {
            if (Double.TryParse(((TextBox)sender).Text, out double tmp))
                b = tmp;
            else
                ((TextBox)sender).Text = b.ToString();
        }

        private void comboBoxFunctions_SelectedIndexChanged(object sender, EventArgs e)
        {
            funcnum = ((ComboBox)sender).SelectedIndex;
        }

        private void trackBar_Scroll(object sender, EventArgs e)
        {
            TrackBar tb = (TrackBar)sender;
            Chart chart;
            tbtochart.TryGetValue(tb, out chart);
            chart.ChartAreas[0].AxisX.ScaleView.ZoomReset();
            chart.ChartAreas[0].AxisY.ScaleView.ZoomReset();
            if (tb.Value == 0)
                return;


            double xfrom = chart.ChartAreas[0].AxisX.ScaleView.ViewMinimum,
                xto = chart.ChartAreas[0].AxisX.ScaleView.ViewMaximum,
                yfrom = chart.ChartAreas[0].AxisY.ScaleView.ViewMinimum,
                yto = chart.ChartAreas[0].AxisY.ScaleView.ViewMaximum;
            double xk = (xto - xfrom) * ((double)tb.TickFrequency / (double)(tb.Maximum * 2+1));
            double yk = (yto - yfrom) * ((double)tb.TickFrequency / (double)(tb.Maximum * 2+1));

            chart.ChartAreas[0].AxisX.ScaleView.Zoom(xfrom + tb.Value * xk, xto - tb.Value * xk);
            chart.ChartAreas[0].AxisY.ScaleView.Zoom(yfrom + tb.Value * yk, yto - tb.Value * yk);

        }

        private void textBoxX_TextChanged(object sender, EventArgs e)
        {
            if(!Double.TryParse(((TextBox)sender).Text, out double x))
                return;
            if (s == null)
                return;


            richTextBoxOutputs.Text = "";
            richTextBoxOutputs.Text += "f(x)    == " + s.OriginFunc(x).ToString("G13") + "\n";
            richTextBoxOutputs.Text += "S(A,x)  == " + s.Eval(x).ToString("G13") + "\n";
            richTextBoxOutputs.Text += "f`(x)   == " + s.OriginDer(x).ToString("G13") + "\n";
            richTextBoxOutputs.Text += "S`(A,x) == " + s.EvalDer(x).ToString("G13") + "\n";

        }

        private void chart_MouseMove(object sender, MouseEventArgs e)
        {
            Chart chart = (Chart)sender;
            try
            {
                double x = chart.ChartAreas[0].AxisX.PixelPositionToValue((double)chart.PointToClient(MousePosition).X);
                double y = chart.ChartAreas[0].AxisY.PixelPositionToValue((double)chart.PointToClient(MousePosition).Y);
                ;
                toolTipPos.SetToolTip(chart, '{' + x.ToString("G5") + ", " + y.ToString("G5") + '}');
            }
            catch
            {
                return;
            }
        }

        private void textBoxNu_Leave(object sender, EventArgs e)
        {
            nu = Double.Parse(((TextBox)sender).Text);
            if (nu == 0.0)
                nu = Double.PositiveInfinity;
        }

        private void buttonEvalAll_Click(object sender, EventArgs e)
        {

            this.Cursor = Cursors.WaitCursor;

            System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            try
            {
                s = new Spline(funcnum, linknum, a, b, nu);
            }
            catch (Exception exc)
            {
                MessageBox.Show(exc.Message, "Помилка", MessageBoxButtons.OK, MessageBoxIcon.Error);
                this.Cursor = Cursors.Default;
                return;
            }
            sw.Stop();
            // richTextBoxOutputs.Text = sw.Elapsed.ToString() + '\n';
            this.Cursor = Cursors.AppStarting;

            dataGridViewParams.Columns["Column4"].Visible = s.param_count == 5;
            dataGridViewParams.Rows.Clear();
            dataGridViewParams.Rows.Add(s.link_count);
            for (int i = 0; i < s.link_count; i++)
            {
                dataGridViewParams.Rows[i].Cells["Num"].Value = (i+1).ToString();
                dataGridViewParams.Rows[i].Cells["Borders"].Value = "[" + s.X[i].ToString("0.00") + ", " + s.X[i + 1].ToString("0.00") + "]";
                dataGridViewParams.Rows[i].Cells["MaxError"].Value = s.MaxError(s.X[i], s.X[i + 1]).ToString("0.00000000E+00");
                //MaxError
                for (int j = 0; j < s.param_count; j++)
                {
                    dataGridViewParams.Rows[i].Cells[j + 3].Value = s.A[i, j].ToString("0.0000E+00");
                }
            }

            chart1.Series[0].Points.Clear();
            chart1.Series[1].Points.Clear();
            chart2.Series[0].Points.Clear();
            chart2.Series[1].Points.Clear();
            chart3.Series[0].Points.Clear();
            chart3.Series[1].Points.Clear();
            chart4.Series[0].Points.Clear();
            chart4.Series[1].Points.Clear();


            //Color[] cs = { Color.Red, Color.Green, Color.Blue };

            for (int k = 0; k < s.link_count; k++) { 
                double step = (s.X[k+1] - s.X[k]) / (2000);
                for (int i = 0; i <= 2000; i++)
                {
                    double x = s.X[k] + i * step;
                    double y1 = s.OriginFunc(x);
                    double y2 = s.Eval(x);
                    double dy1 = s.OriginDer(x);
                    double dy2 = s.EvalDer(x);
                    chart1.Series[0].Points.AddXY(x, y1);
                    chart1.Series[1].Points.AddXY(x, y2);
                    chart2.Series[0].Points.AddXY(x, y1 - y2);
                    chart3.Series[0].Points.AddXY(x, dy1);
                    chart3.Series[1].Points.AddXY(x, dy2);
                    chart4.Series[0].Points.AddXY(x, dy1 - dy2);
                }
            }
            
            foreach (double x in s.X)
            {
                double y1 = chart2.Series[0].Points.FindByValue(x, "X").YValues[0];
                chart2.Series[1].Points.AddXY(x,y1);
                double y2 = chart4.Series[0].Points.FindByValue(x, "X").YValues[0];
                chart4.Series[1].Points.AddXY(x, y2);
            }

            this.Cursor = Cursors.Default;
        }

        private void comboBoxLinks_SelectedIndexChanged(object sender, EventArgs e)
        {
            switch (((ComboBox)sender).SelectedIndex + 1)
            {
                case 1:
                    linknum = (int)linktype.powexp4;
                    pictureBoxLink.Image = HermiteSpline.Properties.Resources.ES4;
                    break;
                case 2:
                    linknum = (int)linktype.powexp5;
                    pictureBoxLink.Image = HermiteSpline.Properties.Resources.ES5;
                    break;
                case 3:
                    linknum = (int)linktype.poly4;
                    pictureBoxLink.Image = HermiteSpline.Properties.Resources.Pl4;
                    break;
                case 4:
                    linknum = (int)linktype.poly5;
                    pictureBoxLink.Image = HermiteSpline.Properties.Resources.Pl5;
                    break;
                case 5:
                    linknum = (int)linktype.exppow5;
                    pictureBoxLink.Image = HermiteSpline.Properties.Resources.EP5;
                    break;
                case 6:
                    linknum = (int)linktype.pow1exp2;
                    pictureBoxLink.Image = HermiteSpline.Properties.Resources.W12;
                    break;
                case 7:
                    linknum = (int)linktype.pow2exp2;
                    pictureBoxLink.Image = HermiteSpline.Properties.Resources.W22;
                    break;
            }
        }
    }
}
