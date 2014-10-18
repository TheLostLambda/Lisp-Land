;; See Licence for legal information

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

;;(defun dump-ci () Work in progress, dump all cell info in a formatted form.
;;  (dotimes (i (length *cells*))
;;  (format t "~%Cell: ~a~%Position: ~a~%Nucleic Acids:"))
