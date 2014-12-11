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

(defun GeneR (celli var)
  (cond ((equal var 'chlorop)
    (let ((Chlorop 0))
      (incf Chlorop (G0 celli var)) Chlorop))))
