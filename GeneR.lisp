(defun G0 (celli var)
  (let ((Gval (parse-gene (nth 0 (fetch-value :DNA celli *cells*))))
        (Chlorop nil))
    (cond ((btwn Gval 0 9)
           (setf Chlorop 1))
          ((btwn Gval 10 19)
           (setf Chlorop 2))
          ((btwn Gval 20 29)
           (setf Chlorop 3))
          (t 
           (setf Chlorop 1))) ;;This is the "default value". AKA if the gene is not expressed.
    (cond ((equal var 'chlorop)
            Chlorop))))
           
(defun G1 (celli var)
  (let ((Gval (parse-gene (nth 1 (fetch-value :DNA celli *cells*))))
        (Perm nil))
    (cond ((btwn Gval 0 9)
           (setf Perm 2))
          ((btwn Gval 10 19)
           (setf Perm 3))
          ((btwn Gval 20 29)
           (setf Perm 4))
          (t 
           (setf Perm 1))) ;;This is the "default value". AKA if the gene is not expressed.
    (cond ((equal var 'perm)
            Perm))))

(defun G2 (celli var)
  (let ((Gval (parse-gene (nth 2 (fetch-value :DNA celli *cells*))))
        (Perm nil)
        (Rigd nil))
    (cond ((btwn Gval 0 7)
           (setf Perm 3)
           (setf Rigd 3))
          ((btwn Gval 8 14)
           (setf Perm 2)
           (setf Rigd 1))
          ((btwn Gval 15 21)
           (setf Perm 1)
           (setf Rigd 1))
          ((btwn Gval 22 28)
           (setf Perm 0)
           (setf Rigd 0))
          (t 			  
           (setf Perm 0) ;;This is the "default value". AKA if the gene is not expressed.
           (setf Rigd 0))) 
    (cond ((equal var 'perm)
            Perm)
          ((equal var 'rigd)
            Rigd))))

(defun GeneR (celli var)
  (cond ((equal var 'chlorop)
          (let ((Chlorop 0))
            (incf Chlorop (G0 celli var)) Chlorop))
        ((equal var 'perm)
          (let ((Perm 0))
            (incf Perm (G1 celli var))
            (incf Perm (G2 celli var)) Perm))
        ((equal var 'rigd)
          (let ((Rigd 0))
            (incf Rigd (G2 celli var)) Rigd))))
