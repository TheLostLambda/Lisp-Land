;; See Licence for legal information

(defparameter *cells* '())

(defclass std-cell ()
  ((POS
    :initarg :POS
    :accessor POS)
  (ATP
    :initarg :ATP
    :accessor ATP)
  (NA 
    :initarg :NA
    :accessor NA)
  (AA
    :initarg :AA
    :accessor AA)
  (FA
    :initarg :FA
    :accessor FA)
  (G
    :initarg :G
    :accessor G)
  (DNA
    :initarg :DNA
    :accessor DNA)))

(defun new-cell (POS ATP NA AA FA G DNA)
  (push (make-instance 'std-cell :POS POS :ATP ATP :NA NA :AA AA 
                                :FA FA :G G :DNA DNA) *cells*))

(defmacro fetch-value (accessor index)
  `(,accessor (nth ,index (reverse *cells*))))

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
