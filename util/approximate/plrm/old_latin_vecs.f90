    REAL(KIND=RT), DIMENSION(SIZE(COLUMN_VECTORS,2)) :: ORDERING
    REAL(KIND=RT) :: NUMBER, TEMP

    ! For each component, create a random ordering of points to assign
    !  them all to bins (to guarantee decent well-spacedness).
    DO I = 1, SIZE(COLUMN_VECTORS, 1)
       ! Create a random ordering that is 0-indexed.
       FORALL (J = 0 : SIZE(COLUMN_VECTORS, 2)-1)
          ORDERING(J) = REAL(J, KIND=RT)
       END FORALL
       DO J = SIZE(COLUMN_VECTORS, 2), 1, -1
          CALL RANDOM_NUMBER(NUMBER)
          NUMBER = 1.0 + NUMBER * REAL(J, KIND=RT)
          IF (NUMBER .LT. J) THEN
             TEMP = ORDERING(INT(NUMBER))
             ORDERING(INT(NUMBER)) = ORDERING(J)
             ORDERING(J) = TEMP
          END IF
       END DO
       ! Use the ordering as the "cell offset" for generated random numbers.
       COLUMN_VECTORS(I,:) = COLUMN_VECTORS(I,:) + ORDERING(:)
       ! Repeat the above for the temporary vectors.
       FORALL (J = 0 : SIZE(COLUMN_VECTORS, 2)-1)
          ORDERING(J) = REAL(J, KIND=RT)
       END FORALL
       DO J = SIZE(COLUMN_VECTORS, 2), 1, -1
          CALL RANDOM_NUMBER(NUMBER)
          NUMBER = 1.0 + NUMBER * REAL(J, KIND=RT)
          IF (NUMBER .LT. J) THEN
             TEMP = ORDERING(INT(NUMBER))
             ORDERING(INT(NUMBER)) = ORDERING(J)
             ORDERING(J) = TEMP
          END IF
       END DO
       TEMP_VECS(I,:) = TEMP_VECS(I,:) + ORDERING(:)
    END DO
    COLUMN_VECTORS(:,:) = COLUMN_VECTORS(:,:) / &
         REAL(SIZE(COLUMN_VECTORS, 2), KIND=RT)
