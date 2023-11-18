using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

namespace HermiteSpline
{
    enum linktype
    {
        powexp4  = 1,
        powexp5  = 2,
        poly4    = 3,
        poly5    = 4,
        exppow5  = 5,
        pow1exp2 = 6,
        pow2exp2 = 7,
    }

    internal class Spline
    {
        internal readonly int funcnum, linknum, param_count, link_count;
        internal readonly double[,] A;
        internal readonly double[] X;


        struct herm_params
        {
            internal int type;
            internal int param_count; //к-сть параметрів у одній ланці
            internal int link_count;  //кількість ланок
            internal IntPtr A;        //параметри сплайна (link_conut*param_count елементів)
            internal IntPtr X;		  //точки наближення (link_count+1 елементів)
            internal IntPtr A128;     //точки наближення (80, 96 або 128 біт)
            internal IntPtr X128;     //параметри сплайна (80, 96 або 128 біт)
        };
        herm_params hp;

        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern herm_params _HermGenNu(Byte funcnum, Byte linknum, double a, double b, double nu);
        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern herm_params _HermGenR(Byte funcnum, Byte linknum, double a, double b, Int32 r);

        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern void _free(herm_params hp);

        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern double _HermiteSpline(herm_params hp, double x, Byte derivative = 0);

        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern double _Func(Byte funcnum, double x, Byte derivative = 0);

        [DllImport("HermiteLib.dll", CallingConvention = CallingConvention.Cdecl)]
        static extern double _MaxError(herm_params hp, Byte funcnum, double from, double to);


        private Spline(herm_params hp)
        {
            this.hp = hp;
            if (hp.X == IntPtr.Zero || hp.A == IntPtr.Zero)
            {
                throw new Exception("Couldn't create spline");
            }
            this.linknum = hp.type;
            this.link_count = this.hp.link_count;
            this.param_count = this.hp.param_count;
            this.A = new double[this.link_count, this.param_count];
            this.X = new double[this.link_count + 1];

            Marshal.Copy(hp.X, X, 0, link_count + 1);
            double[] tmp = new double[link_count * param_count];
            Marshal.Copy((IntPtr)hp.A, tmp, 0, link_count * param_count);
            for (int i = 0; i < link_count; i++)
            {
                for (int j = 0; j < param_count; j++)
                {
                    A[i, j] = tmp[i * param_count + j];
                }
            }
        }

        public Spline(int funcnum, int linknum, double a, double b, double nu) : this(_HermGenNu((Byte)funcnum, (Byte)linknum, a, b, nu))
        {
            /*
            herm_params hp;
            try
            {
                hp = ;
            }
            catch(Exception ex)
            {
                System.Windows.Forms.MessageBox.Show(ex.Message);
                return;
            }
            new Spline(hp);
            */
            this.funcnum = funcnum;
        }
        public Spline(int funcnum, int linknum, double a, double b, int r) : this(_HermGenR((Byte)funcnum, (Byte)linknum, a, b, r))
        {
            /*
            herm_params hp;
            try
            {
                hp = _HermGenR((Byte)funcnum, (Byte)linknum, a, b, r);
            }
            catch (Exception ex)
            {
                System.Windows.Forms.MessageBox.Show(ex.Message);
                return;
            }
            new Spline(hp);
            */
            this.funcnum = funcnum;
        }

        ~Spline()
        {
            _free(this.hp);
            this.hp = default;  //не впевнений навіщо я це роблю
        }

        public double Eval(double x)
        {
            return _HermiteSpline(hp, x);
        }
        public double OriginFunc(double x)
        {
            return _Func((Byte)funcnum, x);
        }
        public double EvalDer(double x)
        {
            return _HermiteSpline(hp, x, 1);
        }
        public double OriginDer(double x)
        {
            return _Func((Byte)funcnum, x, 1);
        }
        public double MaxError(double from, double to)
        {
            return _MaxError(hp, (Byte)funcnum, from, to);
        }

    }
}
