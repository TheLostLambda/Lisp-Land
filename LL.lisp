(defparameter *cells* '())

(defun new-cell (POS ATP NA AA FA G DNA)
  (push (list :POS POS :ATP ATP :NA NA 
         :AA AA :FA FA :G G :DNA DNA) *cells*))

(defmacro fetch-value (accessor index)
  `(getf (nth ,index (reverse *cells*)) ,accessor))

(defun parse-gene (gene)
  (let ((gene-val 0))
    (dotimes (i (length gene))
      (setf gene-val (+ (* (expt 2 i) (nth i (reverse gene))) gene-val))) gene-val))
      
(defun parse-DNA (DNA-seq)
  (let ((DNA-val nil))
    (dotimes (i (length DNA-seq))
      (push (parse-gene (nth i DNA-seq)) DNA-val)) (reverse DNA-val)))

(defun dump-ci ()
  (dolist (entry *cells*)
    (format t "~%~{~a:~10t~a~%~}~%" entry)))

(defun save-ci (filename)
  (with-open-file (out filename
                   :direction :output
                   :if-exists :supersede)
    (with-standard-io-syntax
      (print *cells* out))))

(defun load-ci (filename)
  (with-open-file (in filename)
    (with-standard-io-syntax
      (setf *cells* (read in)))))
      
(defun close-sim ()
  (format t "Saving simulation state and exiting...")
  (save-ci "Celldata.db")
  (quit))
  
(defun init-sim ()
  (format t "Loading simulation state and starting...~%")
  (load-ci "Celldata.db")
  (format t "...done."))
