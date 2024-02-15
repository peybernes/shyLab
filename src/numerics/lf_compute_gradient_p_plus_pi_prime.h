      RealType p_plus_pi_prime_ooo  = p_plus_pi_prime[cell_ooo];  
      RealType p_plus_pi_prime_m1m1 = p_plus_pi_prime[cell_m1m1];
      RealType p_plus_pi_prime_p1p1 = p_plus_pi_prime[cell_p1p1];
      RealType p_plus_pi_prime_m1p1 = p_plus_pi_prime[cell_m1p1];
      RealType p_plus_pi_prime_p1m1 = p_plus_pi_prime[cell_p1m1];
      RealType p_plus_pi_prime_om1  = p_plus_pi_prime[cell_om1];
      RealType p_plus_pi_prime_op1  = p_plus_pi_prime[cell_op1];
      RealType p_plus_pi_prime_m1o  = p_plus_pi_prime[cell_m1o];
      RealType p_plus_pi_prime_p1o  = p_plus_pi_prime[cell_p1o];



      gradx = twelth * (p_plus_pi_prime_p1p1 - p_plus_pi_prime_m1p1) +
	      third  * (p_plus_pi_prime_p1o  - p_plus_pi_prime_m1o ) +
	      twelth * (p_plus_pi_prime_p1m1 - p_plus_pi_prime_m1m1);

      grady = twelth * (p_plus_pi_prime_p1p1 - p_plus_pi_prime_p1m1) +
	      third  * (p_plus_pi_prime_op1  - p_plus_pi_prime_om1 ) +
	      twelth * (p_plus_pi_prime_m1p1 - p_plus_pi_prime_m1m1);


      grad_tmp = 0.5 * grady + 0.5*gradx;
      if(grad_tmp > 0)
	{
	  zpp     = std::max(p_plus_pi_prime_ooo,p_plus_pi_prime_p1o);
	  zpp_tmp = std::max(p_plus_pi_prime_op1,p_plus_pi_prime_p1p1);
	  zpp = std::max(zpp,zpp_tmp);
	  //zsup=1.;
	  zpp = 0.5 * (zpp - p_plus_pi_prime_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zpp     = std::min(p_plus_pi_prime_ooo,p_plus_pi_prime_p1o);
	  zpp_tmp = std::min(p_plus_pi_prime_op1,p_plus_pi_prime_p1p1);
	  zpp=std::min(zpp,zpp_tmp);
      //zsup=0;
	  zpp = 0.5 * (zpp - p_plus_pi_prime_ooo) / (grad_tmp);
	}
      grad_tmp = -0.5 * grady + 0.5 * gradx;
      if(grad_tmp > 0)
	{
	  zpm     = std::max(p_plus_pi_prime_ooo,p_plus_pi_prime_p1o);
	  zpm_tmp =std::max(p_plus_pi_prime_om1,p_plus_pi_prime_p1m1);
	  zpm=std::max(zpm,zpm_tmp);
	  //zsup=1.;
	  zpm = 0.5 * (zpm - p_plus_pi_prime_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zpm     = std::min(p_plus_pi_prime_ooo,p_plus_pi_prime_p1o);
	  zpm_tmp = std::min(p_plus_pi_prime_om1,p_plus_pi_prime_p1m1);
	  zpm = std::min(zpm,zpm_tmp);
      //zsup=0;
	  zpm = 0.5 * (zpm - p_plus_pi_prime_ooo) / (grad_tmp);
	}

      grad_tmp = 0.5 * grady - 0.5 * gradx;
      if(grad_tmp > 0)
	{
	  zmp     = std::max(p_plus_pi_prime_ooo,p_plus_pi_prime_m1o);
	  zmp_tmp = std::max(p_plus_pi_prime_op1,p_plus_pi_prime_m1p1);
	  zmp = std::max(zmp,zmp_tmp);
	  //zsup=1.;
	  zmp = 0.5 * (zmp - p_plus_pi_prime_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zmp     = std::min(p_plus_pi_prime_ooo,p_plus_pi_prime_m1o);
	  zmp_tmp = std::min(p_plus_pi_prime_op1,p_plus_pi_prime_m1p1);
	  zmp = std::min(zmp,zmp_tmp);
      //zsup=0;
	  zmp = 0.5 * (zmp - p_plus_pi_prime_ooo) / (grad_tmp);
	}

      grad_tmp = -0.5 * grady - 0.5 * gradx;
      if(grad_tmp > 0)
	{
	  zmm     = std::max(p_plus_pi_prime_ooo,p_plus_pi_prime_m1o);
	  zmm_tmp = std::max(p_plus_pi_prime_om1,p_plus_pi_prime_m1m1);
	  zmm = std::max(zmm,zmm_tmp);
	  //zsup=1.;
	  zmm = 0.5 * (zmm - p_plus_pi_prime_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zmm     = std::min(p_plus_pi_prime_ooo,p_plus_pi_prime_m1o);
	  zmm_tmp = std::min(p_plus_pi_prime_om1,p_plus_pi_prime_m1m1);
	  zmm = std::min(zmm,zmm_tmp);
      //zsup=0;
	  zmm = 0.5 * (zmm - p_plus_pi_prime_ooo) / (grad_tmp);
	}

      tmp = std::min(zpp,zmm);
      phi = std::min(zmp,zpm);
      phi = std::min(phi,tmp);
      phi = std::min(phi,2.);
