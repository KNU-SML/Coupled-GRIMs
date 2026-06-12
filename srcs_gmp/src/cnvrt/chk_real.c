#include <stdio.h>
#include <math.h>
#include <fp_class.h>

int		chk_real_(data,new,n)
float	*data,*new;
int		*n;
{
	register	i,res;
	float		*fp;

	for (fp=data; fp< data+ *n; fp++)
	{
		res=fp_classf(*fp);
		if (res==FP_SNAN       || res==FP_QNAN ||
			res==FP_POS_DENORM || res==FP_NEG_DENORM)
		{
			/*
			fprintf(stderr,"%d-th data is not valid:%d\n",i,res);
			*/
			*fp= *new;
		}
	}
	return(0);
}
