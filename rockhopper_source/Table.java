/*
 * Copyright 2014 Brian Tjaden
 *
 * This file is part of Rockhopper.
 *
 * Rockhopper is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * Rockhopper is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * (in the file gpl.txt) along with Rockhopper.  
 * If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Hashtable implementation. Uses two parallel arrays. 
 * Uses open addressing for collision resolution.
 */

import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLongArray;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicInteger;

public class Table {

    public int capacity;
    private AtomicInteger size = new AtomicInteger();  // Size of table
    private double loadFactor = 0.90;
    private AtomicLongArray keys;
    private AtomicIntegerArray values;
    private long prime = 16777619;
    private Long offset = new Long("2166136261");
    private int hashPower = (int)Math.pow(2, Assembler.k) - 1;

    public Table() {
	if (Assembler.CAPACITY_POWER == 25) capacity = 30000001;
	else if (Assembler.CAPACITY_POWER >= 31) capacity = (int)(Math.pow(2, 31) - 7);  // Max int
	else if (Assembler.CAPACITY_POWER % 4 == 0) capacity = (int)(Math.pow(2, Assembler.CAPACITY_POWER) + 1);
	else capacity = (int)(Math.pow(2, Assembler.CAPACITY_POWER) - 1);
	keys = new AtomicLongArray(capacity);
	values = new AtomicIntegerArray(capacity);
    }

    /**
     * If "key" is not in table, add it with "value" of 1.
     * If "key" is in table, increment its "value".
     */
    public void add(Long key) {
	int index = hash(key);
	while (values.get(index) != 0) {
	    if (keys.get(index) == key.longValue()) {
		values.incrementAndGet(index);
		return;
	    }
	    index = (index+1) % capacity;
	}
	keys.set(index, key.longValue());
	values.incrementAndGet(index);
	size.incrementAndGet();
    }

    /**
     * Used for adding strand ambiguous reads.
     * If "key" is not in table, add it with "value" of 1.
     * If "key" is in table, increment its "value".
     */
    public void add(Long key, Long key_RC) {
	int index = hash(key);
	while (values.get(index) != 0) {
	    if (keys.get(index) == key.longValue()) {
		values.incrementAndGet(index);
		return;
	    }
	    index = (index+1) % capacity;
	}
	int index_RC = hash(key_RC);
	while (values.get(index_RC) != 0) {
	    if (keys.get(index_RC) == key_RC.longValue()) {
		values.incrementAndGet(index_RC);
		return;
	    }
	    index_RC = (index_RC+1) % capacity;
	}
	keys.set(index, key.longValue());
	values.incrementAndGet(index);
	size.incrementAndGet();
    }

    public double getLoadFactor() {
	return (size.get()/(double)capacity);
    }

    public boolean exceedsLoadFactor() {
	return ((size.get()/(double)capacity) >= loadFactor);
    }

    public boolean containsKey(Long key) {
	int index = hash(key);
	while (values.get(index) != 0) {
	    if (keys.get(index) == key.longValue()) return true;
	    index = (index+1) % capacity;
	}
	return false;
    }

    public int size() {
	return size.get();
    }

    public int get(Long key) {
	int index = hash(key);
	while (values.get(index) != 0) {
	    if (keys.get(index) == key.longValue()) return values.get(index);
	    index = (index+1) % capacity;
	}
	return -1;
    }

    public Long getKeyAtIndex(int i) {
	return (new Long(keys.get(i)));
    }

    public int getValueAtIndex(int i) {
	return values.get(i);
    }

    public void remove(Long key) {
	int index = hash(key);
	while (values.get(index) != 0) {
	    if (keys.get(index) == key.longValue()) {
		keys.set(index, 0);
		values.set(index, 0);
		size.decrementAndGet();
		index = (index+1) % capacity;
		while (values.get(index) != 0) {
		    // Re-hash
		    int index2 = hash(new Long(keys.get(index)));
		    while ((index2 != index) && (values.get(index2) != 0)) index2 = (index2+1) % capacity;
		    if ((index2 != index) && (values.get(index2) == 0)) {
			keys.set(index2, keys.get(index));
			values.set(index2, values.get(index));
			keys.set(index, 0);
			values.set(index, 0);
		    }
		    index = (index+1) % capacity;
		}
		return;
	    }
	    index = (index+1) % capacity;
	}
    }

    private int hash(Long key) {

	// FNV-1a
	long hash = offset.longValue();
	long num = key.longValue();
	for (int i=0; i<2; i++) {
	    hash = hash ^ (num&hashPower);
	    hash *= prime;
	    num = num >>> Assembler.k;
	}
	int result = (int)(hash) % capacity;
	if (result < 0) result += capacity;
	return result;
    }

    /*
    private synchronized void rehash() {
	double REHASH_FACTOR = 1.5;
	if (!exceedsLoadFactor()) return;  // Another thread already rehashed
	int capacity2 = (int)(REHASH_FACTOR*capacity + 1.0);
	AtomicLongArray keys2 = new AtomicLongArray(capacity2);
	AtomicIntegerArray values2 = new AtomicIntegerArray(capacity2);
	for (int i=0; i<capacity; i++) {
	    if (values.get(i) != 0) {
		//int i2 = (new Long(keys.get(i))).hashCode() % capacity2;
		//if (i2 < 0) i2 += capacity2;
		int i2 = hash(new Long(keys.get(i)));
		while (values2.get(i2) != 0) i2 = (i2+1) % capacity2;
		keys2.set(i2, keys.get(i));
		values2.set(i2, values.get(i));
	    }
	}
	capacity = capacity2;
	keys = keys2;
	values = values2;
	System.gc();
    }
    */

    public static void main(String[] args) {

	Table t = new Table();
	t.add(new Long(40));
	t.add(new Long(0));
	t.add(new Long(10));
	t.add(new Long(40));
	t.add(new Long(21));
	t.add(new Long(30));
	t.add(new Long(30));
	t.add(new Long(40));
	t.add(new Long(30));
	t.add(new Long(21));
	t.add(new Long(40));
	t.remove(new Long(30));

	for (int i=0; i<t.capacity; i++) System.out.println(t.keys.get(i) + "\t" + t.values.get(i));
    }

}

