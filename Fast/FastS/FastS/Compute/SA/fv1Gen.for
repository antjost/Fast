          if (cutOffExist==1) then
             DistCoef = rop_cutOff(l)
          else
             DistCoef = 1
          end if
          fvv1     = DistCoef*fv1(DistCoef*chi)
