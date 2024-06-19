!! this is needed for the linearized approach of Tamaki et al. 2017
      dist     = max(dlng(l), 1.e-27)
      DistCoef = 1
      if (dist<=cutOffDist) then
         DistCoef = cutOffDist/dist
      end if
