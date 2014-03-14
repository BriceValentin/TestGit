#include <mapping/Mapping_kernels.h>

void TetraMapper3f_apply( const TVecTetra & map_i, const TVecVec4 & map_f, TVecVec3 & out, const TVecVec3 & in )
{
    for (unsigned int i=0;i<out.size();++i)
    {
        out[i] = 
            in[map_i[i][0]] * map_f[i][0] +
            in[map_i[i][1]] * map_f[i][1] +
            in[map_i[i][2]] * map_f[i][2] +
            in[map_i[i][3]] * map_f[i][3];
    }
}
