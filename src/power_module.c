#include "power_module.h"

void compute_timetraces( gridConfiguration *gridCfg, saveData *saveDCfg, int t_int ){

    if ( (t_int % (int)period) == 4 )  {
            //printf( "status: number of oscillation periods: %d (t_int= %d) \n",T_wave,t_int);
            /*printf( "        Poynting-power: z1 = %13.6e, z2 = %13.6e, x1 = %13.6e, x2 = %13.6e, y1 = %13.6e, y2 = %13.6e, (z1+z2+x1+x2+y1+y2)/z1_ref = %13.6e %%\n",
                    power_abs_z1/power_abs_ref, 
                    power_abs_z2/power_abs_ref,
                    power_abs_x1/power_abs_ref, 
                    power_abs_x2/power_abs_ref,
                    power_abs_y1/power_abs_ref, 
                    power_abs_y2/power_abs_ref,
                    (power_abs_x1+power_abs_x2 + power_abs_y1+power_abs_y2 + power_abs_z1+power_abs_z2)/power_abs_ref * 100.
                    );
            timetraces[T_wave][0]   = (double)t_int;
            timetraces[T_wave][1]   = (double)T_wave;
            timetraces[T_wave][2]   = power_abs_z1/power_abs_ref;
            timetraces[T_wave][3]   = power_abs_z2/power_abs_ref;
            timetraces[T_wave][4]   = power_abs_x1/power_abs_ref;
            timetraces[T_wave][5]   = power_abs_x2/power_abs_ref;
            timetraces[T_wave][6]   = power_abs_y1/power_abs_ref;
            timetraces[T_wave][7]   = power_abs_y2/power_abs_ref;*/

        }

}