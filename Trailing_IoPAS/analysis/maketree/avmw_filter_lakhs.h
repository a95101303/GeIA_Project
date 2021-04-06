float *avmw_filter_lakhs(unsigned int pulse_length, unsigned int window, float input_pulse[], float *smooth_pulse)
{
  
  for(unsigned int bin=0;bin<pulse_length;bin++)
    {
      if(bin<window)
	{
	  float tmp_sum = 0;
	  for(unsigned int rr = bin; rr<bin+window; rr++)
	    {
	      tmp_sum += input_pulse[rr]; 
	    } 
	  smooth_pulse[bin] = tmp_sum/(float)(window);
	  
	}
      if(bin>=window &&bin<=(pulse_length-window))
	{
	  float tmp_sum = 0;
	  for(unsigned int rr = bin-(window/2); rr<bin+(window/2); rr++)
	    {
	      tmp_sum += input_pulse[rr]; 

	    } 
	  smooth_pulse[bin] = tmp_sum/(float)window;
	}
      
      if(bin>(pulse_length-window) && bin<=pulse_length)
	{
	  float tmp_sum = 0;
	  for(unsigned int rr = bin-window; rr<bin; rr++)
	    {
	      tmp_sum += input_pulse[rr]; 
	    } 
	  smooth_pulse[bin] = tmp_sum/(float)window;
	  
	}

    }
  return smooth_pulse;

}
