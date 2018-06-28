      IF(param_int(ITYPZONE).eq.3) THEN !2D

         if(param_int(NEQ).eq.5) then
            
            do k = ind_loop_lu(5), ind_loop_lu(6)
               do j = ind_loop_lu(3) - 1, ind_loop_lu(4) + 1

                  lij = indssor(ind_loop_lu(1), j, k)

                  do l = lij-1, lij + ind_loop_lu(2) - ind_loop_lu(1)+1

                     ssortmp(l) = drodm_in(l) + ssor(l)
                     ssortmp(l + v2) = drodm_in(l + v2) + ssor(l + v2)
                     ssortmp(l + v3) = drodm_in(l + v3) + ssor(l + v3)
                     ssortmp(l + v5) = drodm_in(l + v5) + ssor(l + v5)

                     ssor(l) = 0.
                     ssor(l + v2) = 0.
                     ssor(l + v3) = 0.
                     ssor(l + v5) = 0.

                  enddo
               enddo
            enddo

         else !SA

            do k = ind_loop_lu(5), ind_loop_lu(6)
               do j = ind_loop_lu(3) - 1, ind_loop_lu(4) + 1

                  lij = indssor(ind_loop_lu(1), j, k)

                  do l = lij-1, lij + ind_loop_lu(2) - ind_loop_lu(1)+1

                     ssortmp(l) = drodm_in(l) + ssor(l)
                     ssortmp(l + v2) = drodm_in(l + v2) + ssor(l + v2)
                     ssortmp(l + v3) = drodm_in(l + v3) + ssor(l + v3)
                     ssortmp(l + v5) = drodm_in(l + v5) + ssor(l + v5)
                     ssortmp(l + v6) = drodm_in(l + v6) + ssor(l + v6)

                     ssor(l) = 0.
                     ssor(l + v2) = 0.
                     ssor(l + v3) = 0.
                     ssor(l + v5) = 0.
                     ssor(l + v6) = 0.

                  enddo
               enddo
            enddo

         endif

      else !3D

         if(param_int(NEQ).eq.5) then

            do k = ind_loop_lu(5) - 1, ind_loop_lu(6) + 1
               do j = ind_loop_lu(3) - 1, ind_loop_lu(4) + 1

                  lij = indssor(ind_loop_lu(1), j, k)

                  do l = lij-1, lij + ind_loop_lu(2) - ind_loop_lu(1)+1

                     ssortmp(l) = drodm_in(l) + ssor(l)
                     ssortmp(l + v2) = drodm_in(l + v2) + ssor(l + v2)
                     ssortmp(l + v3) = drodm_in(l + v3) + ssor(l + v3)
                     ssortmp(l + v4) = drodm_in(l + v4) + ssor(l + v4)
                     ssortmp(l + v5) = drodm_in(l + v5) + ssor(l + v5)

                     ssor(l) = 0.
                     ssor(l + v2) = 0.
                     ssor(l + v3) = 0.
                     ssor(l + v4) = 0.
                     ssor(l + v5) = 0.

                  enddo
               enddo
            enddo

         else !SA

            do k = ind_loop_lu(5) - 1, ind_loop_lu(6) + 1
               do j = ind_loop_lu(3) - 1, ind_loop_lu(4) + 1

                  lij = indssor(ind_loop_lu(1), j, k)

                  do l = lij-1, lij + ind_loop_lu(2) - ind_loop_lu(1)+1

                     ssortmp(l) = drodm_in(l) + ssor(l)
                     ssortmp(l + v2) = drodm_in(l + v2) + ssor(l + v2)
                     ssortmp(l + v3) = drodm_in(l + v3) + ssor(l + v3)
                     ssortmp(l + v4) = drodm_in(l + v4) + ssor(l + v4)
                     ssortmp(l + v5) = drodm_in(l + v5) + ssor(l + v5)
                     ssortmp(l + v6) = drodm_in(l + v6) + ssor(l + v6)

                     ssor(l) = 0.
                     ssor(l + v2) = 0.
                     ssor(l + v3) = 0.
                     ssor(l + v4) = 0.
                     ssor(l + v5) = 0.
                     ssor(l + v6) = 0.

                  enddo
               enddo
            enddo

         endif
      endif
