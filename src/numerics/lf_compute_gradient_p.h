    RealType pressure_ooo  = pressure[cell_ooo];  
    RealType pressure_m1m1 = pressure[cell_m1m1];
    RealType pressure_p1p1 = pressure[cell_p1p1];
    RealType pressure_m1p1 = pressure[cell_m1p1];
    RealType pressure_p1m1 = pressure[cell_p1m1];
    RealType pressure_om1  = pressure[cell_om1];
    RealType pressure_op1  = pressure[cell_op1];
    RealType pressure_m1o  = pressure[cell_m1o];
    RealType pressure_p1o  = pressure[cell_p1o];

    tmp_gamma_r = 0;
    tmp_pi_r    = 0;
    tmp_gamma_l = 0;
    tmp_pi_l    = 0;
    tmp_gamma_t = 0;
    tmp_pi_t    = 0;
    tmp_gamma_b = 0;
    tmp_pi_b    = 0;
    rho_e_gradx_left [cell_ooo] = 0;
    rho_e_grady_bot  [cell_ooo] = 0;
    rho_e_gradx_right[cell_ooo] = 0;
    rho_e_grady_top  [cell_ooo] = 0;
    double tmpgradx_rho = 0;
    double tmpgrady_rho = 0;
    zpp = 1.;
    zpm = 1.;
    zmp = 1.;
    zmm = 1.;
      //
      gradx = twelth * (pressure_p1p1 - pressure_m1p1) +
	      third  * (pressure_p1o  - pressure_m1o ) +
	      twelth * (pressure_p1m1 - pressure_m1m1) ;

      grady = twelth * (pressure_p1p1 - pressure_p1m1) +
	      third  * (pressure_op1  - pressure_om1 ) +
	      twelth * (pressure_m1p1 - pressure_m1m1) ;


      grad_tmp = 0.5 * grady + 0.5 * gradx;
      if(grad_tmp > 0)
	{
	  zpp     = std::max(pressure_ooo,pressure_p1o);
	  zpp_tmp = std::max(pressure_op1,pressure_p1p1);
	  zpp     = std::max(zpp,zpp_tmp);
	  //zsup=1.;
	  zpp = (zpp - pressure_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zpp     = std::min(pressure_ooo,pressure_p1o);
	  zpp_tmp = std::min(pressure_op1,pressure_p1p1);
	  zpp     = std::min(zpp,zpp_tmp);
      //zsup=0;
	  zpp = (zpp - pressure_ooo) / (grad_tmp);
	}
      grad_tmp = -0.5 * grady + 0.5 * gradx;
      if(grad_tmp > 0)
	{
	  zpm     = std::max(pressure_ooo,pressure_p1o);
	  zpm_tmp = std::max(pressure_om1,pressure_p1m1);
	  zpm     =std::max(zpm,zpm_tmp);
	  //zsup=1.;
	  zpm = (zpm - pressure_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zpm     = std::min(pressure_ooo,pressure_p1o);
	  zpm_tmp = std::min(pressure_om1,pressure_p1m1);
	  zpm     = std::min(zpm,zpm_tmp);
      //zsup=0;
	  zpm = (zpm - pressure_ooo) / (grad_tmp);
	}

      grad_tmp = 0.5 * grady - 0.5 * gradx;
      if(grad_tmp> 0 )
	{
	  zmp     = std::max(pressure_ooo,pressure_m1o);
	  zmp_tmp = std::max(pressure_op1,pressure_m1p1);
	  zmp     = std::max(zmp,zmp_tmp);
	  //zsup=1.;
	  zmp = (zmp - pressure_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zmp     = std::min(pressure_ooo,pressure_m1o);
	  zmp_tmp =std::min(pressure_op1,pressure_m1p1);
	  zmp     = std::min(zmp,zmp_tmp);
      //zsup=0;
	  zmp = (zmp - pressure_ooo) / (grad_tmp);
	}

      grad_tmp = -0.5 * grady - 0.5 * gradx;
      if(grad_tmp > 0)
	{
	  zmm     = std::max(pressure_ooo,pressure_m1o);
	  zmm_tmp = std::max(pressure_om1,pressure_m1m1);
	  zmm     = std::max(zmm,zmm_tmp);
	  //zsup=1.;
	  zmm = (zmm - pressure_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zmm     = std::min(pressure_ooo,pressure_m1o);
	  zmm_tmp = std::min(pressure_om1,pressure_m1m1);
	  zmm     = std::min(zmm,zmm_tmp);
      //zsup=0;
	  zmm = (zmm - pressure_ooo) / (grad_tmp);
	}
      
      for(int j=0; j < nb_mat; j++){

	RealType alpha_ooo  = in_c_k[j][cell_ooo] ;
	RealType alphak_gradx_left_ooo   = alphak_gradx_left [j][cell_ooo] ;
	RealType alphak_gradx_right_ooo  = alphak_gradx_right[j][cell_ooo] ;
	RealType alphak_grady_top_ooo    = alphak_grady_top  [j][cell_ooo] ;
	RealType alphak_gradx_bot_ooo    = alphak_grady_bot  [j][cell_ooo] ;
	

	tmp_gamma_r = tmp_gamma_r
	  + (alpha_ooo + 0.5 * dx * alphak_gradx_right_ooo)
			       / (gamma_k[j] - 1);

	
	tmp_pi_r = tmp_pi_r
	  + (alpha_ooo + 0.5 * dx * alphak_gradx_right_ooo)
	  * pi_prime_k[j] * gamma_k[j] / (gamma_k[j] - 1);

	tmp_gamma_l = tmp_gamma_l
	  + (alpha_ooo - 0.5 * dx * alphak_gradx_left_ooo)
			       / (gamma_k[j] - 1);

	
	tmp_pi_l = tmp_pi_l
	  + (alpha_ooo - 0.5 * dx * alphak_gradx_left_ooo)
	  * pi_prime_k[j] * gamma_k[j] / (gamma_k[j] - 1);

	tmp_gamma_t = tmp_gamma_t
	  + (alpha_ooo + 0.5 * dy * alphak_grady_top_ooo)
			       / (gamma_k[j] - 1);

	
	tmp_pi_t = tmp_pi_t
	  + (alpha_ooo + 0.5 * dy * alphak_grady_top_ooo)
	  * pi_prime_k[j] * gamma_k[j] / (gamma_k[j] - 1);
	
	tmp_gamma_b = tmp_gamma_b
	  + (alpha_ooo - 0.5 * dy * alphak_gradx_bot_ooo)
			       / (gamma_k[j] - 1);

	tmp_pi_b = tmp_pi_b
	  + (alpha_ooo - 0.5 * dx * alphak_gradx_bot_ooo)
	  * pi_prime_k[j] * gamma_k[j] / (gamma_k[j] - 1);
 
      }
      tmp_gamma_r = 1. / tmp_gamma_r + 1.;
      tmp_gamma_b = 1. / tmp_gamma_b + 1.;
      tmp_gamma_t = 1. / tmp_gamma_t + 1.;
      tmp_gamma_l = 1. / tmp_gamma_l + 1.;
      tmp = std::min(zpp,zmm);
      phi = std::min(zmp,zpm);
      phi = std::min(phi,tmp);
      phi = std::min(0.5*phi,1.0);
      p_grady_top  [cell_ooo]   = phi * grady * h_y;
      p_gradx_right[cell_ooo]   = phi * gradx * h_x;
      p_grady_bot  [cell_ooo]   = phi * grady * h_y;
      p_gradx_left [cell_ooo]   = phi * gradx * h_x;

      if((   pressure[cell_ooo] + 0.5 * dx * p_gradx_right[cell_ooo]) / (tmp_gamma_r - 1) < -tmp_pi_r / tmp_gamma_r
	 || (pressure[cell_ooo] - 0.5 * dx * p_gradx_left [cell_ooo]) / (tmp_gamma_l - 1) < -tmp_pi_l / tmp_gamma_l
	 || (pressure[cell_ooo] + 0.5 * dy * p_grady_top  [cell_ooo]) / (tmp_gamma_t - 1) < -tmp_pi_t / tmp_gamma_t
	 || (pressure[cell_ooo] - 0.5 * dy * p_grady_bot  [cell_ooo]) / (tmp_gamma_b - 1) < -tmp_pi_b / tmp_gamma_b )
	{
	  tmp = std::min(zpp,zmm);
	  phi = std::min(zmp,zpm);
	  phi = std::min(phi,tmp);
	  phi = std::min(phi,2.0);
	  p_grady_top  [cell_ooo] = 2.  * h_y * (pi_prime_mix[cell_ooo] + 0.5 * dy * p_plus_pi_prime_grady[cell_ooo] - tmp_pi_t * (tmp_gamma_t - 1) / tmp_gamma_t);
	  p_grady_bot  [cell_ooo] = -2. * h_y * (pi_prime_mix[cell_ooo] - 0.5 * dy * p_plus_pi_prime_grady[cell_ooo] - tmp_pi_b * (tmp_gamma_b - 1) / tmp_gamma_b);
	  p_gradx_right[cell_ooo] = 2.  * h_x * (pi_prime_mix[cell_ooo] + 0.5 * dx * p_plus_pi_prime_gradx[cell_ooo] - tmp_pi_r * (tmp_gamma_r - 1) / tmp_gamma_r);
	  p_gradx_left [cell_ooo] = -2. * h_x * (pi_prime_mix[cell_ooo] - 0.5 * dx * p_plus_pi_prime_gradx[cell_ooo] - tmp_pi_l * (tmp_gamma_l - 1) / tmp_gamma_l);
	  // grady_p[i]=phi*grady*h_x;
	  //gradx_p[i]=phi*gradx*h_x;
	}

  for(int j=0; j < nb_mat; j++){
      //m_M->setgrady_alpha(k+i,0);
      //m_M->setgradx_alpha(k+i,0);
    RealType alpha_ooo  = in_c_k[j][cell_ooo] ;
    RealType alphak_gradx_left_ooo   = alphak_gradx_left [j][cell_ooo] ;
    RealType alphak_gradx_right_ooo  = alphak_gradx_right[j][cell_ooo] ;
    RealType alphak_grady_top_ooo    = alphak_grady_top  [j][cell_ooo] ;
    RealType alphak_gradx_bot_ooo    = alphak_grady_bot  [j][cell_ooo] ;

    rho_e_gradx_right[cell_ooo] = rho_e_gradx_right[cell_ooo] + (alpha_ooo + 0.5 * dx * alphak_gradx_right_ooo)
      * (pressure[cell_ooo] + 0.5 * dx * p_gradx_right[cell_ooo] + gamma_k[j] * pi_prime_k[j]) / (gamma_k[j] - 1);
    
    rho_e_gradx_left[cell_ooo]  = rho_e_gradx_left [cell_ooo] + (alpha_ooo - 0.5 * dx * alphak_gradx_left_ooo )
      * (pressure[cell_ooo] - 0.5 * dx * p_gradx_left [cell_ooo] + gamma_k[j] * pi_prime_k[j]) / (gamma_k[j] - 1);
    
    rho_e_grady_bot [cell_ooo]  = rho_e_grady_bot  [cell_ooo] + (alpha_ooo - 0.5 * dy * alphak_gradx_bot_ooo  )
      * (pressure[cell_ooo] - 0.5 * dy * p_grady_bot  [cell_ooo] + gamma_k[j] * pi_prime_k[j]) / (gamma_k[j] - 1);

    rho_e_grady_top [cell_ooo]  = rho_e_grady_top  [cell_ooo] + (alpha_ooo + 0.5 * dy * alphak_grady_top_ooo  )
      * (pressure[cell_ooo] + 0.5 * dy * p_grady_top  [cell_ooo] + gamma_k[j] * pi_prime_k[j]) / (gamma_k[j] - 1);
    
  }

    rho_e_gradx_left [cell_ooo]  = -2. * h_x * (rho_e_gradx_left [cell_ooo] - rho_e[cell_ooo]);
    rho_e_grady_bot  [cell_ooo]  = -2. * h_y * (rho_e_grady_bot  [cell_ooo] - rho_e[cell_ooo]);
    rho_e_gradx_right[cell_ooo]  = 2.  * h_x * (rho_e_gradx_right[cell_ooo] - rho_e[cell_ooo]);
    rho_e_grady_top  [cell_ooo]  = 2.  * h_y * (rho_e_grady_top  [cell_ooo] - rho_e[cell_ooo]);      
      // m_M->setgrady_mk(k+i,0);
      //m_M->setgradx_mk(k+i,0);
